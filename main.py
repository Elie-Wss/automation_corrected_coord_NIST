#!/usr/bin/env python
# main.py
import argparse
from io_functions import read_chromato_and_chromato_cube
from processing_functions import read_full_spectra_centroid
from visualization_functions import visualizer

def main():
    parser = argparse.ArgumentParser(description="Process GC GC MS CDF files and produce results.")
    parser.add_argument("samples", nargs="+", help="List of CDF file paths to process")
    args = parser.parse_args()
    
    for sample in args.samples:
        print(f"Processing {sample} ...")
        # Open and process the CDF file
        try:
            chromato, time_rn, chromato_cube, sigma, mass_range = read_chromato_and_chromato_cube(sample)
            # For testing, simply print some info
            print("Chromatogram shape:", chromato.shape)
            print("Time range:", time_rn)
            print("Chromatogram cube shape:", chromato_cube.shape)
            print("Sigma:", sigma)
            print("Mass range:", mass_range)
            
       
        except Exception as e:
            print(f"Failed to process {sample}: {e}")

if __name__ == "__main__":
    main()
