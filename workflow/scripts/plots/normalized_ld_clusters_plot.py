#!/usr/bin/env python

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Normalized LD cluster table (TSV)")
    parser.add_argument("--snp_table", required=True, help="SNP counts per chromosome (TSV)")
    parser.add_argument("--output", required=True, help="Output file prefix (e.g., results/plot_ld_clusters)")
    args = parser.parse_args()

    # Load data
    df_ld = pd.read_csv(args.input, sep='\t')
    df_snp = pd.read_csv(args.snp_table, sep='\t')

    # Extract chromosome numbers
    df_snp['chr'] = df_snp['Chromosome'].str.extract(r'\.chr(\d+)$')[0]
    df_ld['chr'] = df_ld['chr'].str.replace('chr', '').astype(str)

    # Merge and sort
    df = pd.merge(df_ld, df_snp[['chr', 'PASS_SNP_Count']], on='chr', how='left')
    df['chr_num'] = df['chr'].astype(int)
    df = df.sort_values('chr_num')

    # Melt for grouped barplot
    plot_df = pd.melt(
        df,
        id_vars='chr',
        value_vars=['clusters_per_Mb', 'PASS_SNP_Count'],
        var_name='Metric',
        value_name='Value'
    )

    # Normalize SNP counts for visualization (optional: comment out if not needed)
    plot_df['Value'] = plot_df.apply(
        lambda row: row['Value'] / 1e6 if row['Metric'] == 'PASS_SNP_Count' else row['Value'],
        axis=1
    )
    plot_df['Metric'] = plot_df['Metric'].replace({
        'clusters_per_Mb': 'Clusters per Mb',
        'PASS_SNP_Count': 'SNPs (millions)'
    })

    # Plot
    sns.set(style="whitegrid")
    plt.figure(figsize=(12, 6))
    ax = sns.barplot(data=plot_df, x='chr', y='Value', hue='Metric', palette='muted')

    ax.set_title("LD Cluster Density and SNP Count by Chromosome")
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("Value")
    plt.xticks(rotation=45)
    plt.legend(title='Metric', loc='upper right')

    # Save
    plt.tight_layout()
    plt.savefig(f"{args.output}.png", dpi=300)
    plt.savefig(f"{args.output}.pdf")
    plt.savefig(f"{args.output}.svg")

if __name__ == "__main__":
    main()
