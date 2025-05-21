import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D

# -----------------------------
# COLORBLIND-FRIENDLY SETTINGS
# -----------------------------
sns.set(style="whitegrid", context="talk", palette="colorblind")
cb_palette = sns.color_palette("colorblind")
female_color = cb_palette[0]
male_color = cb_palette[3]
diff_color = cb_palette[2]
region_color = cb_palette[6]
arrow_color = cb_palette[0]

# -----------------------------
# ARGPARSE
# -----------------------------
parser = argparse.ArgumentParser(description="Combine heterozygosity plots with gene annotations")
parser.add_argument("--vcf_tab", required=True, help="VCF tabular file with genotypes")
parser.add_argument("--gff", required=True, help="GFF3 annotation file")
parser.add_argument("--seqid", required=True, help="Chromosome/scaffold name")
parser.add_argument("--region_start", type=int, required=True, help="Start of region")
parser.add_argument("--region_end", type=int, required=True, help="End of region")
parser.add_argument("--top_snps", required=False, help="TSV file with top SNPs to highlight")
parser.add_argument("--out_png", required=True)
parser.add_argument("--out_pdf", required=True)
parser.add_argument("--out_svg", required=True)
args = parser.parse_args()

# -----------------------------
# LOAD DATA
# -----------------------------
df = pd.read_csv(args.vcf_tab, sep="\t")
col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
gff = pd.read_csv(args.gff, sep="\t", comment="#", names=col_names)

# -----------------------------
# FUNCTIONS
# -----------------------------
def compute_heterozygosity(gt):
    if not isinstance(gt, str) or gt in ("./.", ".|."):
        return np.nan
    alleles = gt.replace("|", "/").split("/")
    if len(alleles) != 2:
        return np.nan
    return 1 if alleles[0] != alleles[1] else 0

def infer_sex(name):
    name = name.upper()
    if name.startswith("F") or name.startswith("RF"):
        return "F"
    elif name.startswith("M") or name.startswith("RM"):
        return "M"
    else:
        return "U"

def extract_gene_id(attr):
    for field in attr.split(";"):
        if "Name=" in field:
            return field.split("Name=")[-1]
        elif "gene=" in field:
            return field.split("gene=")[-1]
        elif "ID=" in field:
            return field.split("ID=")[-1]
    return "unknown"

# -----------------------------
# PROCESSING
# -----------------------------
gt_columns = [col for col in df.columns if col.endswith(".GT")]
sex_map = {col: infer_sex(col) for col in gt_columns}
female_inds = [col for col in gt_columns if sex_map[col] == "F"]
male_inds = [col for col in gt_columns if sex_map[col] == "M"]
positions = df["POS"]

# Raw heterozygosity
female_het = df[female_inds].applymap(compute_heterozygosity).mean(axis=1)
male_het = df[male_inds].applymap(compute_heterozygosity).mean(axis=1)

# Smoothed heterozygosity
window_size = 50
female_het_smooth = female_het.rolling(window=window_size, min_periods=1).mean()
male_het_smooth = male_het.rolling(window=window_size, min_periods=1).mean()
diff_het = female_het_smooth - male_het_smooth

# -----------------------------
# GENE ANNOTATIONS
# -----------------------------
region = gff[
    (gff["seqid"] == args.seqid) &
    (gff["type"] == "gene") &
    (gff["start"] <= args.region_end) &
    (gff["end"] >= args.region_start)
].copy()
region["gene_id"] = region["attributes"].apply(extract_gene_id)

label_dict = {
    "gene-BLAG_LOCUS17194": "HAO1",
    "gene-BLAG_LOCUS17195": "FLT1"
}
region["gene_label"] = region["gene_id"].apply(lambda gid: label_dict.get(gid, gid))

region_sorted = region.sort_values("start").reset_index(drop=True)
region_sorted["track"] = 0
track_ends = []

for i, row in region_sorted.iterrows():
    placed = False
    for track_idx, end_pos in enumerate(track_ends):
        if row["start"] > end_pos + 1000:
            region_sorted.at[i, "track"] = track_idx
            track_ends[track_idx] = row["end"]
            placed = True
            break
    if not placed:
        region_sorted.at[i, "track"] = len(track_ends)
        track_ends.append(row["end"])

# -----------------------------
# PLOTTING WITH CUSTOM LAYOUT
# -----------------------------
fig = plt.figure(figsize=(14, 13))
gs = gridspec.GridSpec(3, 1, height_ratios=[2, 3, 1])

ax3 = fig.add_subplot(gs[0])  # Raw heterozygosity (top)
ax1 = fig.add_subplot(gs[1], sharex=ax3)  # Smoothed het
ax2 = fig.add_subplot(gs[2], sharex=ax3)  # Gene annotation

# X range
x_min = args.region_start - 5000
x_max = args.region_end + 25000
cutoff_x = args.region_end + 5000

mask = positions <= cutoff_x

# --- Raw Heterozygosity ---
ax3.scatter(positions[mask], female_het[mask], color=female_color, alpha=0.6, s=8)
ax3.scatter(positions[mask], male_het[mask], color=male_color, alpha=0.6, s=8)
ax3.axvspan(args.region_start, args.region_end, color=region_color, alpha=0.1)
ax3.set_ylabel("Raw Het")
ax3.grid(True)
ax3.get_legend().remove() if ax3.get_legend() else None
ax3.tick_params(labelbottom=False)

# --- Smoothed Heterozygosity ---
ax1.plot(positions, female_het_smooth, color=female_color)
ax1.plot(positions, male_het_smooth, color=male_color)
ax1.plot(positions, diff_het, color=diff_color, linestyle="--")
ax1.axvspan(args.region_start, args.region_end, color=region_color, alpha=0.2)
ax1.set_ylabel("Smoothed Het")
ax1.grid(True)
ax1.get_legend().remove() if ax1.get_legend() else None
ax1.tick_params(labelbottom=False)

# --- Gene Annotation Track ---
ax2.axvspan(args.region_start, args.region_end, color=region_color, alpha=0.2)
for _, row in region_sorted.iterrows():
    y = row["track"] * 0.5
    start, end, strand = row["start"], row["end"], row["strand"]
    name = row["gene_label"][:20]
    arrow = mpatches.FancyArrow(
        start if strand == "+" else end,
        y,
        end - start if strand == "+" else start - end,
        0,
        width=0.1,
        length_includes_head=True,
        head_width=0.2,
        head_length=500,
        color=arrow_color,
        alpha=0.9
    )
    ax2.add_patch(arrow)
    ax2.text((start + end) / 2, y + 0.2, name, ha="center", va="bottom", fontsize=7)

ax2.set_ylim(-0.5, max(region_sorted["track"]) * 0.6 + 1)
ax2.set_yticks([])
ax2.set_xlabel("Base Pair Position (bp)")
ax2.grid(True)
# Leave legend alone so Top SNP legend appears if needed

# Apply shared x-limits
for ax in (ax1, ax2, ax3):
    ax.set_xlim(x_min, x_max)

# -----------------------------
# SNP TAGGING (Vertical red bars in ax2)
# -----------------------------
if args.top_snps:
    top_snps = pd.read_csv(args.top_snps, sep="\t")

    # Normalize CHROM values to match args.seqid
    top_snps["CHROM"] = top_snps["CHROM"].replace({"chr4": args.seqid})

    # Filter SNPs to region and chromosome
    top_snps = top_snps[
        (top_snps["CHROM"] == args.seqid) &
        (top_snps["POS"] >= x_min) &
        (top_snps["POS"] <= x_max)
    ]

    # Plot red vertical bars on gene annotation panel
    for _, row in top_snps.iterrows():
        ax2.axvline(x=row["POS"], color="red", linewidth=2.5, linestyle="-", alpha=0.9)

    if not top_snps.empty:
        legend_elements = [
            Line2D([0], [0], color='red', lw=3, label='Top SNP')
        ]
        ax2.legend(handles=legend_elements, loc='upper right')

# Final layout & export
plt.tight_layout(h_pad=1.5)
plt.savefig(args.out_png)
plt.savefig(args.out_pdf)
plt.savefig(args.out_svg)
