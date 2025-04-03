#!/usr/bin/env python
import argparse
import os
from io_functions import read_chromato_and_chromato_cube
from processing_functions import read_full_spectra_centroid, full_spectra_to_chromato_cube
from visualization_functions import visualizer

def process_sample(sample_path, mod_time=1.25, pre_process=True):
    print(f"Processing sample: {sample_path}")
    # Open and process the CDF file
    chromato, time_rn, chromato_cube, sigma, mass_range = read_chromato_and_chromato_cube(
        sample_path, mod_time=mod_time, pre_process=pre_process
    )
    print("Chromatogram shape:", chromato.shape)
    print("Time range:", time_rn)
    print("Chromatogram cube shape:", chromato_cube.shape)
    print("Sigma:", sigma)
    print("Mass range:", mass_range)
    
    # (Future processing steps will be added here)
    
    # Optional: visualize the result (commented out for now)
    # visualizer((chromato, time_rn, chromato_cube), title=f"Sample: {os.path.basename(sample_path)}")
    
    return

def process_directory(directory, mod_time=1.25, pre_process=True):
    # Loop over all files in the directory
    for filename in os.listdir(directory):
        # Only process files with a .cdf extension
        if filename.lower().endswith('.cdf'):
            sample_path = os.path.join(directory, filename)
            try:
                process_sample(sample_path, mod_time=mod_time, pre_process=pre_process)
            except Exception as e:
                print(f"Error processing {sample_path}: {e}")

def main():
    parser = argparse.ArgumentParser(description="Process GC GC MS CDF files to produce results.")
    parser.add_argument("directory", help="Directory containing CDF files to process")
    parser.add_argument("--mod_time", type=float, default=1.25, help="Modulation time for processing")
    parser.add_argument("--no_preprocess", action="store_true", help="Disable baseline pre-processing")
    args = parser.parse_args()
    
    process_directory(args.directory, mod_time=args.mod_time, pre_process=not args.no_preprocess)

if __name__ == "__main__":
    main()
