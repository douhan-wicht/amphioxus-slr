import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# -----------------------------
# STYLE: Colorblind Palette
# -----------------------------
sns.set(style="whitegrid", context="talk", palette="colorblind")
colors = sns.color_palette("colorblind", n_colors=10)

# -----------------------------
# ARGPARSE
# -----------------------------
parser = argparse.ArgumentParser(description="Manhattan plot using Sex_g values across genome")
parser.add_argument("--input", required=True, help="Input CSV with region metrics (e.g., candidates.csv)")
parser.add_argument("--out_png", required=True, help="Output PNG path")
parser.add_argument("--out_pdf", required=True, help="Output PDF path")
parser.add_argument("--out_svg", required=True, help="Output SVG path")
args = parser.parse_args()

# -----------------------------
# LOAD DATA
# -----------------------------
df = pd.read_csv(args.input)

# Extract coordinates
df[['region_start', 'region_end']] = df['region'].str.extract(r':(\d+)-(\d+)').astype(int)
df['region_length'] = df['region_end'] - df['region_start']
df['midpoint'] = df['region_start'] + df['region_length'] // 2

# Chromosome numeric sort key
df['chr_num'] = df['chr'].str.replace("chr", "").astype(int)
df = df.sort_values(['chr_num', 'midpoint'])

# Cumulative position calculation
chromosomes = df['chr'].unique()
chrom_offsets = {}
offset = 0
cumulative_pos = []

for chrom in sorted(chromosomes, key=lambda x: int(x.replace("chr", ""))):
    chr_df = df[df['chr'] == chrom]
    chrom_offsets[chrom] = offset
    cumulative_pos.extend(chr_df['midpoint'] + offset)
    offset += chr_df['region_end'].max() + 1e6  # 1Mb padding

df['cumulative_pos'] = cumulative_pos

# -----------------------------
# PLOTTING
# -----------------------------
fig, ax = plt.subplots(figsize=(14, 6))

# Alternate colors per chromosome
for i, chrom in enumerate(sorted(df['chr'].unique(), key=lambda x: int(x.replace("chr", "")))):
    chrom_data = df[df['chr'] == chrom]
    ax.scatter(
        chrom_data['cumulative_pos'],
        chrom_data['Sex_g'],
        color=colors[i % len(colors)],
        s=20,
        label=chrom if i % 2 == 0 else "",  # label every other chr
        alpha=0.7
    )

# Axis tweaks
ax.set_xlabel("Genomic Position (by chromosome)")
ax.set_ylabel("Sex_g (misgrouping proportion)")
ax.set_title("Manhattan Plot: Sex Association (Sex_g)", weight="bold")
ax.grid(True, axis="y")

# Chromosome labels
chr_labels = []
chr_ticks = []
for chrom in sorted(df['chr'].unique(), key=lambda x: int(x.replace("chr", ""))):
    chr_df = df[df['chr'] == chrom]
    center = chr_df['cumulative_pos'].median()
    chr_labels.append(chrom)
    chr_ticks.append(center)

ax.set_xticks(chr_ticks)
ax.set_xticklabels(chr_labels, rotation=45)

# Invert y-axis for low Sex_g at top
ax.invert_yaxis()

plt.tight_layout()
plt.savefig(args.out_png)
plt.savefig(args.out_pdf)
plt.savefig(args.out_svg)
