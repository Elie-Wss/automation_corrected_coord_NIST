# annotation.py
import os
import numpy as np
import pandas as pd
import logging
from skimage.feature import peak_local_max
from visualization_functions import chromato_to_matrix, matrix_to_chromato
from processing_functions import read_spectrum_from_chromato_cube
import pyms_nist_search
import pyms

# -------------------------------
# Docker / NIST Engine Setup Helper
# -------------------------------
def ensure_pyms_nist_container_clean():
    import docker
    client = docker.from_env()
    try:
        container = client.containers.get("pyms-nist-server")
        if container.status == "exited":
            print("🧹 Removing existing stopped container: pyms-nist-server", flush=True)
            container.remove(force=True)
    except docker.errors.NotFound:
        pass

# -------------------------------
# Coordinate Correction
# -------------------------------
def correct_coordinates(chromato_cube, chromato, time_rn, mod_time, ref_row, chromato_shape):
    name = ref_row['Mol']
    mass = int(ref_row['mass'])
    RT1 = float(ref_row['RT1']) * 60
    RT2 = float(ref_row['RT2'])
    theoretical_coord = chromato_to_matrix(np.array([[RT1 / 60, RT2]]), time_rn, mod_time, chromato_shape)[0]

    rt1_window = 0.5
    rt2_window = 0.2
    window = chromato_to_matrix(
        np.array([
            [RT1 / 60 - rt1_window, RT2 - rt2_window],
            [RT1 / 60 + rt1_window, RT2 + rt2_window]
        ]),
        time_rn,
        mod_time,
        chromato.shape
    )

    x_min = int(max(float(window[0][0]), 0))
    y_min = int(max(float(window[0][1]), 0))
    x_max = int(min(float(window[1][0]), float(chromato.shape[0])))
    y_max = int(min(float(window[1][1]), float(chromato.shape[1])))

    mass_index = int(mass) - 39
    submatrix = chromato_cube[mass_index, x_min:x_max, y_min:y_max]
    local_max_coords = peak_local_max(submatrix, num_peaks=1)

    if len(local_max_coords) > 0:
        peak_x, peak_y = local_max_coords[0]
        global_x = x_min + peak_x
        global_y = y_min + peak_y
        coord = [global_x, global_y]
        rt_corr = matrix_to_chromato(np.array([coord]), time_rn, mod_time, chromato_shape)[0]
        RT1_corr = rt_corr[0] * 60
        RT2_corr = rt_corr[1]
        return RT1_corr, RT2_corr, coord
    else:
        print(f"⚠️ No local max found near {ref_row['Mol']}, using theoretical RT", flush=True)
        return RT1, RT2, theoretical_coord

# -------------------------------
# Spectrum Extraction & NIST Search
# -------------------------------
def run_nist_annotation(chromato_cube, time_rn, chromato_shape, mass, corrected_coord, mass_range, search_engine, n_hits=20):
    int_values = read_spectrum_from_chromato_cube(corrected_coord, chromato_cube=chromato_cube)
    range_min, range_max = mass_range
    mass_values = np.linspace(range_min, range_max, range_max - range_min + 1).astype(int)
    mass_spectrum = pyms.Spectrum.MassSpectrum(mass_values, int_values)
    return search_engine.full_search_with_ref_data(mass_spectrum, n_hits=n_hits)

# -------------------------------
# Hit Selection – PATCHED VERSION
# -------------------------------
def select_nist_hit(results, variations):
    selected_hit = None
    match_flag = 0
    for hit in results:
        if hit and hit[0]:
            df = pd.DataFrame.from_dict(hit[0]).T
            hit_data = df.iloc[1].to_dict() if df.shape[0] >= 2 else df.iloc[0].to_dict()
            nist_name = hit_data.get('NIST_Name', "").lower()
            if any(variation.lower() in nist_name for variation in variations):
                return hit_data, 1
    # fallback to top hit
    if results and results[0] and results[0][0]:
        df = pd.DataFrame.from_dict(results[0][0]).T
        selected_hit = df.iloc[1].to_dict() if df.shape[0] >= 2 else df.iloc[0].to_dict()
    return selected_hit, 0

# -------------------------------
# QC Metric
# -------------------------------
def compute_intensity_noise_ratio(int_values, sigma):
    return np.max(int_values) / sigma if sigma > 0 else np.nan

# -------------------------------
# Main Annotation Function
# -------------------------------
def annotate_sample_nist(sample_name, chromato, time_rn, chromato_cube, sigma, mass_range,
                         ref_peaks_csv, output_dir, search_engine, canonical_dict):
    ref_peaks = pd.read_csv(ref_peaks_csv)
    fame_peaks = ref_peaks.iloc[9:].reset_index(drop=True)
    results = []
    chromato_shape = chromato.shape

    for idx, row in fame_peaks.iterrows():
        try:
            name = row['Mol']
            theoretical_mass = int(row['mass'])
            RT1_theoretical = float(row['RT1']) * 60
            RT2_theoretical = float(row['RT2'])
            print(f"🔍 Processing reference peak {idx+1}/{len(fame_peaks)}: {name}", flush=True)

            RT1_corr, RT2_corr, corrected_coord = correct_coordinates(
                chromato_cube, chromato, time_rn, 1.25, row, chromato_shape
            )

            nist_results = run_nist_annotation(
                chromato_cube, time_rn, chromato_shape,
                theoretical_mass, corrected_coord, mass_range,
                search_engine, n_hits=50
            )

            variations = canonical_dict.get(name, [])
            selected_hit, match_flag = select_nist_hit(nist_results, variations)
            int_values = read_spectrum_from_chromato_cube(corrected_coord, chromato_cube=chromato_cube)
            intensity_noise_ratio = compute_intensity_noise_ratio(int_values, sigma)

            result_row = {
                'Sample': sample_name,
                'Mol': name,
                'mass': theoretical_mass,
                'RT1_theoretical': RT1_theoretical / 60,
                'RT2_theoretical': RT2_theoretical,
                'RT1_corrected': RT1_corr / 60,
                'RT2_corrected': RT2_corr,
                #'match_flag': match_flag,
                'intensity_noise_ratio': intensity_noise_ratio,
            }

            if selected_hit:
                for key, value in selected_hit.items():
                    result_row[f"NIST_{key}"] = value
            else:
                result_row['NIST_Name'] = None

            results.append(result_row)
        except Exception as e:
            print(f"⚠️ Error processing reference peak {idx} ({row['Mol']}): {e}", flush=True)

    if results:
        return pd.DataFrame(results)
    else:
        print("⚠️ No results to write for sample:", sample_name, flush=True)
        return None
    
def compute_final_match_flag(row):
    """
    Compute the final match flag for a given row.
    Uses the stored NIST_0 field and the molecule name.
    """
    mol = row["Mol"]
    nist_name = str(row.get("NIST_0", "")).lower().strip()
    variations = canonical_dict.get(mol, [])
    return 1 if any(variation.lower().strip() in nist_name for variation in variations) else 0   

