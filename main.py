import os
import argparse
from io_functions import read_chromato_and_chromato_cube
from visualization_functions import visualizer
from pdf_generation import generate_ref_peaks_pdf  # if you saved it in pdf_generation.py

def main():
    parser = argparse.ArgumentParser(description="Process GC GC MS CDF files and produce results.")
    parser.add_argument("directory", help="Directory containing CDF files to process")
    parser.add_argument("--mod_time", type=float, default=1.25, help="Modulation time for processing")
    parser.add_argument("--no_preprocess", action="store_true", help="Disable baseline pre-processing")
    args = parser.parse_args()
    
    for filename in os.listdir(args.directory):
        if filename.lower().endswith('.cdf'):
            sample_path = os.path.join(args.directory, filename)
            print(f"Processing {sample_path} ...")
            try:
                chromato, time_rn, chromato_cube, sigma, mass_range = read_chromato_and_chromato_cube(
                    sample_path, mod_time=args.mod_time, pre_process=not args.no_preprocess
                )
                print("Chromatogram shape:", chromato.shape)
                print("Time range:", time_rn)
                print("Chromatogram cube shape:", chromato_cube.shape)
                print("Sigma:", sigma)
                print("Mass range:", mass_range)
                
                # Extract sample name (without extension) to use in the PDF filename
                sample_name = os.path.splitext(filename)[0]
                generate_ref_peaks_pdf(chromato_cube, time_rn, sample_name, mod_time=args.mod_time)
                
            except Exception as e:
                print(f"Failed to process {sample_path}: {e}")

if __name__ == "__main__":
    main()

