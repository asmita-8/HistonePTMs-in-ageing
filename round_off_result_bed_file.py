import os
import pandas as pd

## Function to round values to the nearest multiple of 200
def round_to_nearest_200(value):
    return int(round(value / 200.0) * 200)

## Function to process a single BED file
def process_bed_file(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            parts = line.strip().split()
            if len(parts) >= 3:
                chrom = parts[0]
                start = round_to_nearest_200(int(parts[1]))
                end = round_to_nearest_200(int(parts[2]))
                outfile.write(f"{chrom}\t{start}\t{end}\n")

## Function to process all BED files in the directory
def process_all_bed_files(base_dir):
    for gene_dir in os.listdir(base_dir):
        gene_path = os.path.join(base_dir, gene_dir)
        if os.path.isdir(gene_path):  ## Ensure it's a directory
            input_bed_file = os.path.join(gene_path, f"{gene_dir}.bed.result.bed")
            output_bed_file = os.path.join(gene_path, f"{gene_dir}_rounded_result.bed")
            
            if os.path.exists(input_bed_file):
                process_bed_file(input_bed_file, output_bed_file)
                print(f"Processed: {input_bed_file} -> {output_bed_file}")
            else:
                print(f"Skipping {gene_dir}, no BED file found.")

## Base directory containing gene subdirectories - parent directory having one directory per gene
base_directory = "/home/asmita/my_dir1/nucpos_200/NucPosSimulator_linux64/nucpos_3yr_200/h3k4me3heartRV34yr_nucpos_results"
process_all_bed_files(base_directory)
print("Processing complete for all genes.")

