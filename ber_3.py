import pandas as pd
import numpy as np

# File paths
file1 = '3_year_final_LV.csv'
file2 = '3yr_final_RV.csv'

# Load data
df1 = pd.read_csv(file1)
df2 = pd.read_csv(file2)

# Extract common gene_ids
common_gene_ids = set(df1['gene_id']).intersection(set(df2['gene_id']))

# All modifications
modifications = ['h3k4me1', 'h3k4me3', 'h3k9me3','h3k27me3', 'h3k36me3']

# Group definitions
groups = {
    'neutral': ['h3k4me1', 'h3k36me3'],
    'activating': ['h3k4me3'],
    'repressive': ['h3k9me3', 'h3k27me3']
}


def calculate_ber_fraction(seq1, seq2):
    total_bits = len(seq1)
    differing_bits = sum(c1 != c2 for c1, c2 in zip(seq1, seq2))
    return differing_bits / total_bits if total_bits > 0 else np.nan


def gene_wise_ber(gene_id):
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

# Process each group separately
for group_name, group_mods in groups.items():
    # Drop any genes that are missing BER values for the required modifications
    group_df = ber_decimal_df[group_mods].dropna()

    # Keep genes where all group modifications have BER > 0.75
    filtered = group_df[(group_df > 0.75).all(axis=1)].copy()

    if filtered.empty:
        print(f"No genes met the >75% BER criteria for group: {group_name}")
        continue

    # Add average column
    filtered["average"] = filtered.mean(axis=1)

    # Sort by average BER descending
    sorted_df = filtered.sort_values(by="average", ascending=False)

    # Save to CSV
    output_path = f"ber_over75_{group_name}_3yr_LVvsRV.csv"
    sorted_df.to_csv(output_path, index=True)
    print(f"Filtered genes for group '{group_name}' saved to: {output_path}")

