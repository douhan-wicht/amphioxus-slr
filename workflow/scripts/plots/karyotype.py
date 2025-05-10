import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib import cm
from matplotlib.colors import Normalize

# -----------------------------
# COLORBLIND-FRIENDLY STYLE
# -----------------------------
sns.set(style="whitegrid", context="talk", palette="colorblind")

# -----------------------------
# ARGPARSE
# -----------------------------
parser = argparse.ArgumentParser(description="Plot karyotype of SLR regions colored by Sex_g")

parser.add_argument("--candidates", required=True, help="Path to candidates.csv")
parser.add_argument("--sex_filter", required=True, help="Path to sex_filter.csv")
parser.add_argument("--out_png", required=True, help="Output PNG path")
parser.add_argument("--out_pdf", required=True, help="Output PDF path")
parser.add_argument("--out_svg", required=True, help="Output SVG path")

args = parser.parse_args()

# -----------------------------
# LOAD DATA
# -----------------------------
candidates = pd.read_csv(args.candidates)
sex_filter = pd.read_csv(args.sex_filter)

# Step 1: Match and merge by common columns
common_cols = candidates.columns.intersection(sex_filter.columns)
df = pd.concat([candidates, sex_filter[common_cols]], ignore_index=True)

# Step 2: Extract start/end from region string
df[['region_start', 'region_end']] = df['region'].str.extract(r':(\d+)-(\d+)').astype(int)

# Step 3: Compute derived columns
df['region_length'] = df['region_end'] - df['region_start']
df['midpoint'] = df['region_start'] + df['region_length'] / 2

# Chromosome lengths dictionary
chr_lengths = {
    '1': 43860960, '2': 38510819, '3': 34610492, '4': 31719604,
    '5': 25701974, '6': 24533633, '7': 24230189, '8': 23752511,
    '9': 23231292, '10': 20381850, '11': 20367708, '12': 19917020,
    '13': 19776172, '14': 19709165, '15': 19381563, '16': 18823661,
    '17': 18214296, '18': 17113871, '19': 15322015
}

# Sort chromosomes numerically
chromosomes = sorted(df['chr'].unique(), key=lambda x: int(x.replace("chr", "")))

# -----------------------------
# PLOTTING (VERTICAL CHROMOSOMES)
# -----------------------------
fig, ax = plt.subplots(figsize=(len(chromosomes) * 1.2, 14))
x_spacing = 2
xticks = []
max_chr_length = max(chr_lengths.values())

# Normalize color scale
norm = Normalize(vmin=df['Sex_g'].min(), vmax=df['Sex_g'].max())
cmap = cm.get_cmap('viridis')  # colorblind-friendly

# Plot chromosomes and regions
for i, chr_name in enumerate(chromosomes):
    x = i * x_spacing
    xticks.append(x)

    chr_data = df[df['chr'] == chr_name]
    chr_num = chr_name.replace("chr", "")
    chr_length = chr_lengths.get(chr_num, chr_data['region_end'].max())

    ax.vlines(x=x, ymin=0, ymax=chr_length, color='black', linewidth=8, alpha=0.5)

    for _, row in chr_data.iterrows():
        color = cmap(norm(row['Sex_g']))
        ax.add_patch(plt.Rectangle(
            (x - 0.5, row['region_start']),
            1.0,
            row['region_length'],
            color=color,
            edgecolor='none',
            alpha=0.95
        ))

# Final layout
ax.set_xticks(xticks)
ax.set_xticklabels(chromosomes, rotation=90, fontsize=10)
ax.set_ylim(0, max_chr_length * 1.05)
ax.set_xlim(-1, max(xticks) + x_spacing)
ax.set_ylabel("Genomic Position (bp)", fontsize=12)
ax.set_title("SLR Regions Colored by Sex_g", fontsize=14, weight='bold')
sns.despine()

# Colorbar
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax, orientation='vertical', pad=0.01)
cbar.set_label("Sex_g (misgrouping proportion)", fontsize=10)

plt.tight_layout()
plt.savefig(args.out_png)
plt.savefig(args.out_pdf)
plt.savefig(args.out_svg)
