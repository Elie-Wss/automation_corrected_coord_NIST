# annotation.py first draft
import os
import re
import logging
import numpy as np
import pandas as pd
import requests
from bs4 import BeautifulSoup
import skimage
from skimage.feature import peak_local_max
from visualization_functions import chromato_to_matrix
from visualization_functions import matrix_to_chromato
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
        pass  # Container doesn't exist
# -------------------------------
# Helper: Coordinate Correction
# -------------------------------
def correct_coordinates(chromato_cube, chromato, time_rn, mod_time, ref_row, chromato_shape):
    name = ref_row['Mol']
    mass = int(ref_row['mass'])
    RT1 = float(ref_row['RT1']) * 60  # seconds
    RT2 = float(ref_row['RT2'])
    u1 = chromato_to_matrix(np.array([[RT1 / 60, RT2]]), time_rn, mod_time, chromato_shape)
    theoretical_coord = u1[0]
    
    rt1_window = 0.5  # minutes
    rt2_window = 0.2  # seconds
    # Define window in matrix coords
    window = chromato_to_matrix(
        np.array([
            [RT1/60 - rt1_window, RT2 - rt2_window],
            [RT1/60 + rt1_window, RT2 + rt2_window]
        ]),
        time_rn,
        1.25,
        chromato.shape
    )

    print("Window:", window, flush=True)
    print("Window shape:", window.shape, flush=True)
    print("window[0] type:", type(window[0]), "value:", window[0], flush=True)
    print("window[1] type:", type(window[1]), "value:", window[1], flush=True)

    # Explicitly extract scalar integer values:
    # Extract scalar bounds from the 2x2 window array explicitly:
    # After computing and printing the window...
    x_min = int(max(float(window[0][0]), 0))
    y_min = int(max(float(window[0][1]), 0))
    x_max = int(min(float(window[1][0]), float(chromato.shape[0])))
    y_max = int(min(float(window[1][1]), float(chromato.shape[1])))
    print("Extracted window values (as ints):", x_min, y_min, x_max, y_max, flush=True)
    print("Final indices:", x_min, y_min, x_max, y_max, flush=True)
    print("Types:", type(x_min), type(y_min), type(x_max), type(y_max), flush=True)

    mass_index = int(mass) - 39
    i_x_min = int(x_min)
    i_y_min = int(y_min)
    i_x_max = int(x_max)
    i_y_max = int(y_max)
    print("Explicit indices:", mass_index, i_x_min, i_x_max, i_y_min, i_y_max, flush=True)

    # Now slice using these explicitly converted integers
    submatrix = chromato_cube[mass_index, i_x_min:i_x_max, i_y_min:i_y_max]
    local_max_coords = peak_local_max(submatrix, num_peaks=1)



        
    if len(local_max_coords) > 0:
        peak_x, peak_y = local_max_coords[0]
        global_x = int(i_x_min) + peak_x
        global_y = int(i_y_min) + peak_y
        coord = [global_x, global_y]
        rt_corr = matrix_to_chromato(np.array([coord]), time_rn, mod_time, chromato_shape)[0]
        RT1_corr = rt_corr[0] * 60  # seconds
        RT2_corr = rt_corr[1]
        return RT1_corr, RT2_corr, coord  # Return the corrected matrix coordinate
    else:
        print(f"⚠️ No local max found near {ref_row['Mol']}, using theoretical RT", flush=True)
        return RT1, RT2, theoretical_coord


# -------------------------------
# Helper: Spectrum Extraction & NIST Search
# -------------------------------
def run_nist_annotation(chromato_cube, time_rn, chromato_shape, mass, corrected_coord, mass_range, search_engine, n_hits=20):
    """
    Extracts a mass spectrum from the chromatogram cube at the given coordinate,
    constructs a mass spectrum object, and runs the NIST search.
    Returns the search results.
    """
    int_values = read_spectrum_from_chromato_cube(corrected_coord, chromato_cube=chromato_cube)
    range_min, range_max = mass_range
    mass_values = np.linspace(range_min, range_max, range_max - range_min + 1).astype(int)
    mass_spectrum = pyms.Spectrum.MassSpectrum(mass_values, int_values)
    
    results = search_engine.full_search_with_ref_data(mass_spectrum, n_hits=n_hits)
    return results

# -------------------------------
# Helper: Hit Selection
# -------------------------------
def select_nist_hit(results, canonical_dict):
    """
    Iterates through all NIST search results for a reference peak and selects the hit
    whose 'NIST_Name' contains one of the acceptable substrings from canonical_dict.
    If none match, falls back to the top hit.
    
    Returns a tuple: (selected_hit, match_flag)
      - selected_hit: A dictionary with NIST hit details.
      - match_flag: 1 if a canonical match was found; 0 otherwise.
    """
    selected_hit = None
    match_flag = 0
    for hit in results:
        if hit and hit[0]:
            df = pd.DataFrame.from_dict(hit[0]).T
            if df.shape[0] >= 2:
                hit_data = df.iloc[1].to_dict()
            else:
                hit_data = df.iloc[0].to_dict()
            nist_name = hit_data.get('NIST_Name', "").lower()
            for canonical, variations in canonical_dict.items():
                if any(variation.lower() in nist_name for variation in variations):
                    selected_hit = hit_data
                    match_flag = 1
                    return selected_hit, match_flag
    if results and results[0] and results[0][0]:
        df = pd.DataFrame.from_dict(results[0][0]).T
        if df.shape[0] >= 2:
            selected_hit = df.iloc[1].to_dict()
        else:
            selected_hit = df.iloc[0].to_dict()
        match_flag = 0
    return selected_hit, match_flag

# -------------------------------
# Helper: QC Metric (Intensity-to-Noise Ratio)
# -------------------------------
def compute_intensity_noise_ratio(int_values, sigma):
    """
    Computes a simple intensity-to-background noise ratio.
    """
    if sigma > 0:
        return np.max(int_values) / sigma
    else:
        return np.nan

# -------------------------------
# Main Annotation Function for a Sample
# -------------------------------
def annotate_sample_nist(sample_name, chromato, time_rn, chromato_cube, sigma, mass_range,
                         ref_peaks_csv, output_dir, search_engine, canonical_dict):
    """
    Processes one sample: for each reference peak in the reference CSV (filtered for FAMEs),
    it performs coordinate correction, spectrum extraction, NIST search, hit selection, and
    computes QC metrics.
    
    The final output is a DataFrame with one row per reference peak for this sample.
    Also writes a sample-specific CSV.
    """
    ref_peaks = pd.read_csv(ref_peaks_csv)
    fame_peaks = ref_peaks.iloc[9:].reset_index(drop=True)
    
    results = []
    chromato_shape = chromato.shape  # (rows, cols)
    for idx, row in fame_peaks.iterrows():
        try:
            name = row['Mol']
            theoretical_mass = int(row['mass'])
            RT1_theoretical = float(row['RT1']) * 60  # seconds
            RT2_theoretical = float(row['RT2'])
            print(f"🔍 Processing reference peak {idx+1}/{len(fame_peaks)}: {name}", flush=True)
            
            RT1_corr, RT2_corr, corrected_coord = correct_coordinates(chromato_cube, chromato, time_rn, 1.25, row, chromato_shape)
            # Use the corrected RT values as the coordinate for spectrum extraction (adjust if needed)
            
            
            nist_results = run_nist_annotation(chromato_cube, time_rn, chromato_shape,
                                   theoretical_mass, corrected_coord, mass_range,
                                   search_engine, n_hits=50)
            selected_hit, match_flag = select_nist_hit(nist_results, canonical_dict)
            int_values = read_spectrum_from_chromato_cube(corrected_coord, chromato_cube=chromato_cube)
            intensity_noise_ratio = compute_intensity_noise_ratio(int_values, sigma)
            
            result_row = {
                'Sample': sample_name,
                'Mol': name,
                'mass': theoretical_mass,
                'RT1_theoretical': RT1_theoretical / 60,  # convert to minutes
                'RT2_theoretical': RT2_theoretical,
                'RT1_corrected': RT1_corr / 60,            # minutes
                'RT2_corrected': RT2_corr,
                'match_flag': match_flag,
                'intensity_noise_ratio': intensity_noise_ratio,
                'NIST_Name': selected_hit.get('NIST_Name') if selected_hit else None,
                'match_factor': selected_hit.get('match_factor') if selected_hit else None
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
        final_df = pd.DataFrame(results)
        return final_df
    else:
        print("⚠️ No results to write for sample:", sample_name, flush=True)
        return None
