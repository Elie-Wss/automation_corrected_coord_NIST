import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from visualization_functions import visualizer

def generate_ref_peaks_pdf(chromato_cube, time_rn, sample_name,
                           ref_peaks_csv="/home/elie/Projects/Work/data_exploration/peak_coordinates.csv",
                           output_dir="/home/elie/Data/test/test_automation_output/pdf_output",
                           mod_time=1.25):
    """
    Generates a PDF of reference chromatogram plots for a given sample.
    
    Parameters:
      chromato_cube: 3D array (mass, rows, cols) of chromatogram data.
      time_rn: Time range tuple.
      sample_name: Name of the sample (used to build the output file name).
      ref_peaks_csv: Path to the CSV file with reference peak coordinates.
      output_dir: Directory where the output PDF will be saved.
      mod_time: Modulation time for processing (passed to the visualizer).
    """
    try:
        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, f"{sample_name}_ref_peaks.pdf")
        print(f"Inside generate_ref_peaks_pdf for sample: {sample_name}", flush=True)
        print(f"Output path: {output_path}", flush=True)
    
        # Load the reference peaks from CSV
        ref_peaks = pd.read_csv(ref_peaks_csv)
        print(f"Loaded {len(ref_peaks)} reference peaks from {ref_peaks_csv}", flush=True)
    
        # Create a multi-page PDF to save each plot
        with PdfPages(output_path) as pdf:
            for i, (_, row) in enumerate(ref_peaks.iterrows(), start=1):
                try:
                    m = int(row['mass'])
                    RT1 = float(row['RT1']) * 60  # convert minutes to seconds
                    RT2 = float(row['RT2'])
                    name = row['Mol']
                    
                    # Adjust mass index (assuming that mass index corresponds to m - 39)
                    m_index = m - 39
                    title = f"{name} | Mass: {m} | RT1: {RT1/60:.2f} min | RT2: {RT2:.3f} s"
                    
                    print(f"[{i}/{len(ref_peaks)}] Plotting: {title}", flush=True)
                    
                    # Call the visualizer function. Adjust the tuple as needed.
                    visualizer(
                        (chromato_cube[m_index, :, :], time_rn),
                        title=title,
                        log_chromato=False,
                        rt1=RT1 / 60,    # convert seconds back to minutes for display
                        rt2=RT2,
                        rt1_window=0.5,
                        rt2_window=0.2,
                        mod_time=mod_time,
                        show=False      # Suppress interactive display
                    )
                    
                    # Save the current figure to the PDF and close it
                    pdf.savefig(plt.gcf())
                    plt.close(plt.gcf())
                    print(f"Saved plot {i}/{len(ref_peaks)} to PDF", flush=True)
                except Exception as e:
                    print(f"Error while plotting reference peak at index {i}: {e}", flush=True)
        print(f"✅ PDF saved to:\n{output_path}", flush=True)
    except Exception as e:
        print(f"Failed in generate_ref_peaks_pdf for sample {sample_name}: {e}", flush=True)

