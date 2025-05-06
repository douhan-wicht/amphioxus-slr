#!/usr/bin/env python3

import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

# -----------------------------
# ARGPARSE
# -----------------------------
parser = argparse.ArgumentParser(description="Plot simplified heatmap of gene expression.")
parser.add_argument("--input", required=True, help="Input TSV file from Bgee with expression values.")
parser.add_argument("--out_png", required=True, help="Output PNG file path.")
parser.add_argument("--out_pdf", required=True, help="Output PDF file path.")
parser.add_argument("--out_svg", required=True, help="Output SVG file path.")
parser.add_argument("--gene", required=True, help="Name of the gene to display in the title.")  # NEW ARGUMENT
args = parser.parse_args()

# -----------------------------
# COLORBLIND-FRIENDLY SETTINGS
# -----------------------------
sns.set(style="whitegrid", context="talk", palette="colorblind")
cb_palette = sns.color_palette("colorblind")
highlight_color = cb_palette[3]  # Red for male-specific

# -----------------------------
# FUNCTION TO SIMPLIFY ANATOMY
# -----------------------------
def simplify_anatomy(name):
    name = name.lower()
    if 'male' in name:
        return 'Testis'  # Renaming Male-specific to Testis
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

# -----------------------------
# LOAD & PROCESS DATA
# -----------------------------
df = pd.read_csv(args.input, sep="\t")
df['category'] = df['anatEntityName'].apply(simplify_anatomy)
df_grouped = df.groupby('category')['value'].mean().reset_index()

# Ordering for readability
category_order = ['Testis', 'Nervous system', 'Musculature', 
                  'Digestive system', 'Developmental stage', 'Other']
df_grouped['category'] = pd.Categorical(df_grouped['category'], categories=category_order, ordered=True)
df_grouped = df_grouped.sort_values('category')

# -----------------------------
# PLOT HEATMAP
# -----------------------------
plt.figure(figsize=(10, 5))  # Adjusted size for better readability
heatmap_data = df_grouped.set_index('category').T

# Create heatmap
ax = sns.heatmap(
    heatmap_data,
    cmap="coolwarm",  # Keeping original color palette
    cbar_kws={'label': 'Expression Level'},
    linewidths=0.5,
    linecolor='gray',
    square=False,  # Adjusting aspect ratio for better layout
    annot=True,  # Show values in cells
    fmt=".1f",  # Formatting the annotation
    annot_kws={"size": 12},  # Adjusting annotation font size
    cbar=True,  # Show colorbar for expression scale
)

# Improve title and axis labels
plt.title(f"{args.gene} expression levels in Branchiostoma lanceolatum", fontsize=16)
plt.xlabel("Tissues", fontsize=14)  # Renaming Tissue Categories to Tissues
plt.ylabel("Expression Level", fontsize=14)

# Rotating x-ticks for better visibility
plt.xticks(rotation=45, ha='right', fontsize=12)
plt.yticks(rotation=0, fontsize=12)

# Remove colorbar label text from the left
cbar = ax.collections[0].colorbar
cbar.set_label('')  # Remove label on the left side of the colorbar

# Adjust layout for better fitting
plt.tight_layout()

# -----------------------------
# EXPORT
# -----------------------------
plt.savefig(args.out_png)
plt.savefig(args.out_pdf)
plt.savefig(args.out_svg)
print(f"Plots saved: {args.out_png}, {args.out_pdf}, {args.out_svg}")
