#!/usr/bin/env python3

import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# -----------------------------
# ARGPARSE
# -----------------------------
parser = argparse.ArgumentParser(description="Plot combined heatmap of two genes.")
parser.add_argument("--input1", required=True, help="Input TSV file for gene 1.")
parser.add_argument("--gene1", required=True, help="Name of gene 1.")
parser.add_argument("--input2", required=True, help="Input TSV file for gene 2.")
parser.add_argument("--gene2", required=True, help="Name of gene 2.")
parser.add_argument("--out_png", required=True, help="Output PNG file path.")
parser.add_argument("--out_pdf", required=True, help="Output PDF file path.")
parser.add_argument("--out_svg", required=True, help="Output SVG file path.")
args = parser.parse_args()

# -----------------------------
# SETTINGS
# -----------------------------
sns.set(style="whitegrid", context="talk", palette="colorblind")

# -----------------------------
# CATEGORY SIMPLIFICATION
# -----------------------------
def simplify_anatomy(name):
    name = name.lower()
    if 'male' in name:
        return 'Testis'
    elif 'nervous' in name:
        return 'Nervous system'
    elif 'musculature' in name or 'muscle' in name:
        return 'Musculature'
    elif 'digestive' in name or 'gut' in name:
        return 'Digestive system'
    elif any(stage in name for stage in ['embryo', 'larva', 'neurula', 'blastula', 'gastrula']):
        return 'Developmental stage'
    else:
        return 'Other'

category_order = ['Testis', 'Nervous system', 'Musculature', 
                  'Digestive system', 'Developmental stage', 'Other']

# -----------------------------
# LOAD AND PROCESS
# -----------------------------
def process_file(filepath, gene_name):
    df = pd.read_csv(filepath, sep="\t")
    df['category'] = df['anatEntityName'].apply(simplify_anatomy)
    df_grouped = df.groupby('category')['value'].mean().reindex(category_order).reset_index()
    df_grouped.columns = ['category', gene_name]
    return df_grouped

df1 = process_file(args.input1, args.gene1)
df2 = process_file(args.input2, args.gene2)

df_merged = pd.merge(df1, df2, on='category')
df_melted = df_merged.melt(id_vars='category', var_name='Gene', value_name='Expression')

# -----------------------------
# PLOT
# -----------------------------
plt.figure(figsize=(10, 5))
pivot_data = df_melted.pivot(index="Gene", columns="category", values="Expression")

ax = sns.heatmap(
    pivot_data,
    cmap="coolwarm",
    annot=True,
    fmt=".1f",
    linewidths=0.5,
    linecolor='gray',
    cbar_kws={'label': 'Expression Score'}
)

plt.title("Expression profiles of candidate genes in Branchiostoma lanceolatum", fontsize=16)
plt.xlabel("Tissues")
plt.ylabel("Genes")
plt.xticks(rotation=45, ha='right')
plt.yticks(rotation=0)
plt.tight_layout()

# -----------------------------
# EXPORT
# -----------------------------
plt.savefig(args.out_png)
plt.savefig(args.out_pdf)
plt.savefig(args.out_svg)
print(f"Saved plots to {args.out_png}, {args.out_pdf}, {args.out_svg}")
