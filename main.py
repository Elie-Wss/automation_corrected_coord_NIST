#!/usr/bin/env python
import os
import argparse
from io_functions import read_chromato_and_chromato_cube
from pdf_generation import generate_ref_peaks_pdf
from annotation import (
    annotate_sample_nist,
    ensure_pyms_nist_container_clean
)
import pyms_nist_search

def process_samples(directory, ref_peaks_csv, annotation_output_dir, pdf_output_dir, mod_time=1.25, pre_process=True):
    # Ensure the output directories exist.
    os.makedirs(annotation_output_dir, exist_ok=True)
    os.makedirs(pdf_output_dir, exist_ok=True)
    
    # Clean up Docker container and initialize NIST search engine.
    ensure_pyms_nist_container_clean()
    search_engine = pyms_nist_search.Engine(
        "/media/elie/Métabolomique/NIST14/MSSEARCH/mainlib/",
        pyms_nist_search.NISTMS_MAIN_LIB,
        "/home/elie/Projects/Work/data_exploration/Python-2DGC-Alignment",
    )
    
    # Define the canonical dictionary for hit selection.
    canonical_dict = {
        "Nonanoic acid, methyl esther C9": ["Nonanoic acid, methyl ester"],
        "Decanoic acid, methyl ester C10": ["Decanoic acid, methyl ester"],
        "Dodecanoic acid, methyl ester C12": ["Dodecanoic acid, methyl ester"],
        "Methyl tetradecanoate C14":["Methyl tetradecanoate"],
        "Hexadecanoic acid, methyl ester C16":["Hexadecanoic acid, methyl ester"],
        "Methyl stearate C18":["Methyl stearate C18"],
        "Eicosanoic acid, methyl ester C20":["Eicosanoic acid, methyl ester"],
        "Docosanoic acid, methyl ester C22":["Docosanoic acid, methyl ester"],
        "Tetracosanoic acid, methyl ester C24":["Tetracosanoic acid, methyl ester"]



    }
    
    all_annotations = []  # List to collect annotation DataFrames from all samples.
    
    for filename in os.listdir(directory):
        if filename.lower().endswith('.cdf'):
            sample_path = os.path.join(directory, filename)
            print(f"Processing sample: {sample_path}")
            try:
                # Read the CDF sample once.
                chromato, time_rn, chromato_cube, sigma, mass_range = read_chromato_and_chromato_cube(
                    sample_path, mod_time=mod_time, pre_process=pre_process
                )
                sample_name = os.path.splitext(filename)[0]
                
                # Run NIST annotation for the sample.
                annotation_df = annotate_sample_nist(
                    sample_name, chromato, time_rn, chromato_cube,
                    sigma, mass_range, ref_peaks_csv, annotation_output_dir,
                    search_engine, canonical_dict
                )
                if annotation_df is not None:
                    all_annotations.append(annotation_df)
                
                # Generate the PDF for the sample.
                generate_ref_peaks_pdf(
                    chromato_cube, time_rn, sample_name,
                    ref_peaks_csv=ref_peaks_csv, output_dir=pdf_output_dir,
                    mod_time=mod_time
                )
            except Exception as e:
                print(f"Failed to process {sample_path}: {e}")
    
    # Aggregate annotation DataFrames from all samples into one master CSV.
    if all_annotations:
        master_df = all_annotations[0].copy()
        for df in all_annotations[1:]:
            master_df = master_df.append(df, ignore_index=True)
        master_csv = os.path.join(annotation_output_dir, "master_nist_annotations.csv")
        master_df.to_csv(master_csv, index=False)
        print(f"✅ Master CSV saved: {master_csv}")
    else:
        print("No annotations generated from any sample.")

def main():
    parser = argparse.ArgumentParser(
        description="Automated NIST Annotation and PDF Generation for GC GC MS Samples"
    )
    parser.add_argument("directory", help="Directory containing CDF sample files")
    parser.add_argument("--ref_peaks_csv", default="/home/elie/Projects/Work/data_exploration/peak_coordinates.csv",
                        help="Path to the reference peaks CSV")
    parser.add_argument("--annotation_output_dir", default="/home/elie/Data/test/test_automation_output/nist_results",
                        help="Directory to save the annotation CSV outputs")
    parser.add_argument("--pdf_output_dir", default="/home/elie/Data/test/test_automation_output/pdf_output",
                        help="Directory to save the generated PDF outputs")
    parser.add_argument("--mod_time", type=float, default=1.25, help="Modulation time for processing")
    parser.add_argument("--no_preprocess", action="store_true", help="Disable baseline pre-processing")
    args = parser.parse_args()
    
    process_samples(
        args.directory,
        args.ref_peaks_csv,
        args.annotation_output_dir,
        args.pdf_output_dir,
        mod_time=args.mod_time,
        pre_process=not args.no_preprocess
    )

if __name__ == "__main__":
    main()


