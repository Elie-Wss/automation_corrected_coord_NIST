import os
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
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, f"{sample_name}_ref_peaks.pdf")
    
    # Load the reference peaks from CSV
    ref_peaks = pd.read_csv(ref_peaks_csv)
    
    # Create a multi-page PDF to save each plot
    with PdfPages(output_path) as pdf:
        for _, row in ref_peaks.iterrows():
            m = int(row['mass'])
            RT1 = float(row['RT1']) * 60  # convert minutes to seconds
            RT2 = float(row['RT2'])
            name = row['Mol']
            
            # Adjust mass index (assuming that mass index corresponds to m - 39)
            m_index = m - 39
            title = f"{name} | Mass: {m} | RT1: {RT1/60:.2f} min | RT2: {RT2:.3f} s"
            
            print(f"📌 Plotting: {title}")
            
            # Call the visualizer function. Adjust the input as needed:
            # Here we assume that visualizer accepts a tuple of (chromato, time_rn, chromato_cube)
            # and other parameters to customize the plot.
            visualizer(
                (chromato_cube[m_index, :, :], time_rn, chromato_cube),
                title=title,
                log_chromato=False,
                rt1=RT1 / 60,    # convert seconds back to minutes for display
                rt2=RT2,
                rt1_window=0.5,
                rt2_window=0.2,
                mod_time=mod_time,
                show=False      # Do not display interactively
            )
            
            # Save the current figure to the PDF and close it
            pdf.savefig(plt.gcf())
            plt.close(plt.gcf())
    
    print(f"✅ PDF saved to:\n{output_path}")
