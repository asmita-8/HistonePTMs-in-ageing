#!/bin/bash

## base directory for the input files - subsetted bed files per gene
input_dir="h3k36me3heartRV34yr_subset_bed_files"

## output parent directory for all gene-specific directories
output_parent_dir="h3k36me3heartRV34yr_nucpos_results"

## Create the parent output directory if it doesn't exist
mkdir -p "$output_parent_dir"

## path to the NucPosSimulator directory - where NucPosSimulator has been downloaded
nucpossimulator_dir="/home/asmita/my_dir1/nucpos_200/NucPosSimulator_linux64/nucpos_3yr_200"  

## Ensure write permissions for the input directory and subdirectories
chmod -R u+w "$input_dir"

## Export variables and define the function for parallel execution
export nucpossimulator_dir output_parent_dir
run_simulator() {
    local bed_file="$1"
    local gene_id
    gene_id=$(basename "$bed_file" .bed)  # Keep _subset in gene_id(just a name)
    
    ## Create a directory for the gene inside the parent output directory
    local gene_dir="$output_parent_dir/$gene_id"
    mkdir -p "$gene_dir"
    
    ## Copy the gene's .bed file to its corresponding directory
    cp "$bed_file" "$gene_dir/"
    
    ## Navigate to the gene directory
    cd "$gene_dir" || exit
    
    ## Run NucPosSimulator in the gene's directory with the params.txt file
    "$nucpossimulator_dir/NucPosSimulator" "$gene_id.bed" "$nucpossimulator_dir/params.txt"
    
    ## Check if the NucPosSimulator ran successfully
    if [[ $? -ne 0 ]]; then
        echo "Error: NucPosSimulator failed for '$gene_id'."
    else
        echo "NucPosSimulator ran successfully for '$gene_id'."
    fi
}

## Export the function for use with GNU Parallel
export -f run_simulator

## Display the number of CPU cores being used
echo "Using 15 CPU cores for parallel execution."

## Find all .bed files in the input directory
find "$input_dir" -type f -name "*_subset.bed" | parallel --jobs 15 --bar run_simulator {}

echo "NucPosSimulator completed for all gene subsets."

