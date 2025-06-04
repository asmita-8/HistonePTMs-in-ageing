import pandas as pd

files = ["discretised_rounded_values_heartRV_h3k4me1_34yr.csv", "discretised_rounded_values_heartRV_h3k4me3_34yr.csv", "discretised_rounded_values_heartRV_h3k27me3_34yr.csv", "discretised_rounded_values_heartRV_h3k36me3_34yr.csv"]
modifications = ["h3k4me1", "h3k4me3","h3k27me3", "h3k36me3"]

## Initialize an empty dataframe
merged_df = None

for file, mod in zip(files, modifications):

    df = pd.read_csv(file)
    
    ## Keep common columns and rename the last column to the corresponding histone mark
    df = df.rename(columns={df.columns[-1]: mod})
    
    if merged_df is None:
        ## First file, initialize the merged dataframe
        merged_df = df
    else:
        ## Merge on common columns
        merged_df = pd.merge(merged_df, df, on=["gene_id", "gene_name", "start_pos", "end_pos", "bin_start", "bin_end"], how="inner")

## Save to a new CSV file
merged_df.to_csv("discretised_rounded_values_34yr_heartRV.csv", index=False)

print("done")
