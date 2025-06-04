import pandas as pd
import numpy as np

# File paths
file1 = '3_year_final_LV.csv'
file2 = '34yr_final_LV.csv'

# Load data
df1 = pd.read_csv(file1)
df2 = pd.read_csv(file2)

# Extract common gene_ids
common_gene_ids = set(df1['gene_id']).intersection(set(df2['gene_id']))

# Histone modifications to consider
modifications = ['h3k4me1', 'h3k4me3','h3k9me3' ,'h3k27me3', 'h3k36me3']

def calculate_ber_fraction(seq1, seq2):
    """Calculate Bit Error Rate as a fraction"""
    total_bits = len(seq1)
    differing_bits = sum(c1 != c2 for c1, c2 in zip(seq1, seq2))
    return differing_bits / total_bits if total_bits > 0 else np.nan

def gene_wise_ber(gene_id):
    """Calculate BER for each modification for a given gene_id."""
    gene_df1 = df1[df1['gene_id'] == gene_id].sort_values(by=['bin_start', 'bin_end'])
    gene_df2 = df2[df2['gene_id'] == gene_id].sort_values(by=['bin_start', 'bin_end'])

    ber_results_decimal = {}
    for mod in modifications:
        seq1 = ''.join(gene_df1[mod].astype(str).tolist())
        seq2 = ''.join(gene_df2[mod].astype(str).tolist())
        ber_results_decimal[mod] = calculate_ber_fraction(seq1, seq2)

    return ber_results_decimal

# Calculate BER for all common genes
ber_decimal_summary = {}
for gene in common_gene_ids:
    decimal_data = gene_wise_ber(gene)
    ber_decimal_summary[gene] = decimal_data

# Convert to DataFrame
ber_decimal_df = pd.DataFrame(ber_decimal_summary).T

# Filter genes where all modifications have BER < 0.05
genes_filtered = ber_decimal_df[ber_decimal_df.max(axis=1) < 0.05]

# Save filtered results
if genes_filtered.empty:
    print("No genes met the filtering criteria.")
else:
    output_file = 'ber_all_mods_less_than_5percent_LV_3vs34.csv'
    genes_filtered.to_csv(output_file, index=True)
    print(f"Filtered genes saved to {output_file}")

