##BER HISTOGRAM FOR ALL GENES
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# File paths
file1 = '/home/ibab/sem_4/final_files/LV_RV/3_year_final_LV.csv'
file2 = '/home/ibab/sem_4/final_files/LV_RV/34yr_final_LV.csv'

# Load data
df1 = pd.read_csv(file1)
df2 = pd.read_csv(file2)

# Extract gene IDs directly from one of the files (assuming both have the same set)
gene_ids = df1['gene_id'].unique()

# Histone modifications to consider
modifications = ['h3k4me1', 'h3k4me3', 'h3k9me3', 'h3k27me3', 'h3k36me3']

def calculate_ber_fraction(seq1, seq2):
    """Calculate Bit Error Rate as a fraction and decimal."""
    total_bits = len(seq1)
    differing_bits = sum(c1 != c2 for c1, c2 in zip(seq1, seq2))
    return f"{differing_bits}/{total_bits}", differing_bits / total_bits if total_bits > 0 else np.nan

def gene_wise_ber(gene_id):
    """Calculate BER as fractions and decimals for each modification for a given gene_id."""
    gene_df1 = df1[df1['gene_id'] == gene_id].sort_values(by=['bin_start', 'bin_end'])
    gene_df2 = df2[df2['gene_id'] == gene_id].sort_values(by=['bin_start', 'bin_end'])

    ber_results_fraction = {}
    ber_results_decimal = {}
    for mod in modifications:
        # Concatenate binary strings for each version
        seq1 = ''.join(gene_df1[mod].astype(str).tolist())
        seq2 = ''.join(gene_df2[mod].astype(str).tolist())
        fraction, decimal = calculate_ber_fraction(seq1, seq2)
        ber_results_fraction[mod] = fraction
        ber_results_decimal[mod] = decimal
    return ber_results_fraction, ber_results_decimal

#Calculate BER for all genes
ber_fraction_summary = {}
ber_decimal_summary = {}

for gene in gene_ids:
    fraction_data, decimal_data = gene_wise_ber(gene)
    ber_fraction_summary[gene] = fraction_data
    ber_decimal_summary[gene] = decimal_data

# DataFrames: Fractions for display, Decimals for visualization & saving
ber_fraction_df = pd.DataFrame(ber_fraction_summary).T
ber_decimal_df = pd.DataFrame(ber_decimal_summary).T

# Save BER results to CSV
ber_fraction_df.to_csv('/home/ibab/sem_4/final_files/LV_RV/ber_all_genes_frac_LV.csv')
ber_decimal_df.to_csv('/home/ibab/sem_4/final_files/LV_RV/ber_all_genes_deci_LV.csv')

#Display summary of BER as fractions
print("Bit Error Rate (BER) as Fractions for all genes:\n")
print(ber_fraction_df.head())

#Histograms per modification
plt.figure(figsize=(15, 10))
bins = 30  # Adjust the number of bins as needed

for idx, mod in enumerate(modifications, 1):
    plt.subplot(2, 3, idx)
    sns.histplot(ber_decimal_df[mod], bins=bins, kde=True, color='blue', edgecolor='black')
    plt.title(f'BER Distribution: {mod}', fontsize=14)
    plt.xlabel('Bit Error Rate (BER)', fontsize=12)
    plt.ylabel('Number of Genes', fontsize=12)
    plt.xlim(0, 1)
    plt.grid(axis='y', linestyle='--', alpha=0.7)

plt.suptitle('Distribution of Bit Error Rate (BER) per Histone Modification', fontsize=16, y=1.02)
plt.tight_layout()
plt.show()

