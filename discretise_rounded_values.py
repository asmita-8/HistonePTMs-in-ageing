import os
import sys
import signal    
import psutil
import pandas as pd
import numpy as np
from pathlib import Path
from joblib import Parallel, delayed
from tqdm import tqdm

## Resource limits
N_JOBS = 50  # 40% of 128 cores
MEMORY_LIMIT = 100 * 1024 * 1024 * 1024  # 100GB
CHUNK_SIZE = 50_000

## Modification configuration
MOD_FOLDER = 'h3k4me3heartRV34yr_nucpos_results'  # Change for different modification
MOD_COL = 'h3k4me3'
BASE_PATH = '/home/asmita/my_dir1/nucpos_200/NucPosSimulator_linux64/nucpos_3yr_200'


def process_chunk(chunk_data):
    """Process a single chunk for the modification."""
    try:
        chunk = chunk_data.copy()
        chunk[MOD_COL] = 0  # Default all values to 0
        
        for gene_id in tqdm(chunk['gene_id'].unique(), leave=False):
            bed_path = Path(BASE_PATH) / MOD_FOLDER / f'{gene_id}_subset' / f'{gene_id}_subset_rounded_result.bed'     
            
            if bed_path.is_file():
                bed_positions = pd.read_csv(
                    bed_path, 
                    sep='\t',
                    usecols=[1, 2],  # Start and end positions
                    header=None,
                    dtype={1: np.int32, 2: np.int32}
                ).values.tolist()
                
                if bed_positions:
                    mask = chunk['gene_id'] == gene_id
                    chunk.loc[mask, MOD_COL] = chunk.loc[mask].apply(
                        lambda row: 1 if [row['bin_start'], row['bin_end']] in bed_positions else 0,
                        axis=1
                    )
        return chunk
    
    except Exception as e:
        print(f"Error processing chunk: {e}")
        return None


def main():
    try:
        ## Read CSV in chunks
        print(f"Processing modification: {MOD_FOLDER}")
        chunks = []
        
        for chunk in tqdm(pd.read_csv('bins_for_all_genes.csv', chunksize=CHUNK_SIZE), 
                         desc="Reading chunks"):
            chunks.append(chunk)
        
        ## Process chunks in parallel
        results = Parallel(n_jobs=N_JOBS, backend='multiprocessing')(
            delayed(process_chunk)(chunk)
            for chunk in tqdm(chunks, desc="Processing chunks")
        )
        
        ## Combine valid results
        valid_results = [r for r in results if r is not None]
        if valid_results:
            final_df = pd.concat(valid_results)
            
            ## Save results
            output_file = f'discretised_rounded_values_heartRV_{MOD_COL}_34yr.csv'
            final_df.to_csv(output_file, index=False)
            print(f"Saved results to {output_file}")
        
    except KeyboardInterrupt:
        print("\nProcessing interrupted by user")
        sys.exit(0)
    except Exception as e:
        print(f"Error in main process: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
