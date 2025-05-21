#!/usr/bin/env python

import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Path to LD cluster table (CSV or TSV)")
    parser.add_argument("--output", required=True, help="Output TSV file for normalized counts")
    args = parser.parse_args()

    # Chromosome lengths (from BraLan3 genome)
    chr_lengths = {
        'chr1': 43860960, 'chr2': 38510819, 'chr3': 34610492, 'chr4': 31719604,
        'chr5': 25701974, 'chr6': 24533633, 'chr7': 24230189, 'chr8': 23752511,
        'chr9': 23231292, 'chr10': 20381850, 'chr11': 20367708, 'chr12': 19917020,
        'chr13': 19776172, 'chr14': 19709165, 'chr15': 19381563, 'chr16': 18823661,
        'chr17': 18214296, 'chr18': 17113871, 'chr19': 15322015
    }

    # Load cluster table
    clusters = pd.read_csv(args.input)
    if 'chr' not in clusters.columns:
        raise ValueError("Input table must contain a 'chr' column")

    # Count clusters per chromosome
    cluster_counts = clusters.groupby('chr').size().reset_index(name='n_clusters')

    # Convert dictionary to DataFrame
    chr_df = pd.DataFrame(list(chr_lengths.items()), columns=['chr', 'length_bp'])

    # Merge and normalize
    merged = pd.merge(cluster_counts, chr_df, on='chr', how='inner')
    merged['clusters_per_Mb'] = merged['n_clusters'] / (merged['length_bp'] / 1e6)

    # Save output
    merged.to_csv(args.output, sep='\t', index=False)

if __name__ == "__main__":
    main()
