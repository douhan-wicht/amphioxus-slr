import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec

# -----------------------------
# ARGPARSE
# -----------------------------
parser = argparse.ArgumentParser(description="Combine heterozygosity plots and gene annotation into one figure")

parser.add_argument("--vcf_tab", required=True, help="VCF tabular file with genotypes")
parser.add_argument("--gff", required=True, help="Path to GFF3 annotation file")
parser.add_argument("--seqid", required=True, help="Chromosome/scaffold name")
parser.add_argument("--region_start", type=int, required=True, help="Start coordinate of region")
parser.add_argument("--region_end", type=int, required=True, help="End coordinate of region")
parser.add_argument("--out_png", required=True, help="Output PNG path")
parser.add_argument("--out_pdf", required=True, help="Output PDF path")
parser.add_argument("--out_svg", required=True, help="Output SVG path")

args = parser.parse_args()

# -----------------------------
# STYLE SETUP
# -----------------------------
sns.set(style="whitegrid", context="talk", palette="colorblind")
palette = sns.color_palette("colorblind")
female_color = palette[0]
male_color = palette[3]
diff_color = palette[2]
region_color = palette[6]
arrow_color = palette[0]

# -----------------------------
# LOAD GENOTYPE DATA
# -----------------------------
df = pd.read_csv(args.vcf_tab, sep="\t")
positions = df["POS"]

def compute_heterozygosity(gt):
    if not isinstance(gt, str) or gt in ("./.", ".|."):
        return np.nan
    alleles = gt.replace("|", "/").split("/")
    return 1 if len(alleles) == 2 and alleles[0] != alleles[1] else 0

def infer_sex(name):
    name = name.upper()
    return "F" if name.startswith("F") or name.startswith("RF") else \
           "M" if name.startswith("M") or name.startswith("RM") else "U"

gt_columns = [col for col in df.columns if col.endswith(".GT")]
sex_map = {col: infer_sex(col) for col in gt_columns}
female_inds = [col for col in gt_columns if sex_map[col] == "F"]
male_inds = [col for col in gt_columns if sex_map[col] == "M"]

female_het = df[female_inds].applymap(compute_heterozygosity).mean(axis=1)
male_het = df[male_inds].applymap(compute_heterozygosity).mean(axis=1)

# Smoothed values
window_size = 50
female_het_smooth = female_het.rolling(window=window_size, min_periods=1).mean()
male_het_smooth = male_het.rolling(window=window_size, min_periods=1).mean()
diff_het = female_het_smooth - male_het_smooth

# -----------------------------
# LOAD GFF & EXTRACT GENES
# -----------------------------
col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
gff = pd.read_csv(args.gff, sep="\t", comment="#", names=col_names)

region = gff[
    (gff["seqid"] == args.seqid) &
    (gff["type"] == "gene") &
    (gff["start"] <= args.region_end) &
    (gff["end"] >= args.region_start)
].copy()

def extract_gene_id(attr):
    for field in attr.split(";"):
        if "Name=" in field:
            return field.split("Name=")[-1]
        elif "gene=" in field:
            return field.split("gene=")[-1]
        elif "ID=" in field:
            return field.split("ID=")[-1]
    return "unknown"

region["gene_id"] = region["attributes"].apply(extract_gene_id)

# Manual label mapping (customize as needed)
label_dict = {
    "BLAG_LOCUS17194": "HAO1",
    "BLAG_LOCUS17195": "FLT1"
}
region["gene_label"] = region["gene_id"].apply(lambda gid: label_dict.get(gid, gid))

# Assign non-overlapping tracks
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
# PLOT FIGURE
# -----------------------------
fig = plt.figure(figsize=(14, 9))
gs = gridspec.GridSpec(nrows=3, ncols=2, width_ratios=[20, 1], height_ratios=[2, 1, 2], hspace=0.3)

# Axes
ax_raw = fig.add_subplot(gs[0, 0])
ax_legend = fig.add_subplot(gs[0, 1])
ax_gene = fig.add_subplot(gs[1, 0])
ax_smooth = fig.add_subplot(gs[2, 0])

# Hide legend axis frame
ax_legend.axis("off")

# X-limits
xlim_raw = (args.region_start - 5000, args.region_end + 5000)
xlim_full = (args.region_start - 5000, args.region_end + 50000)

# -----------------------------
# Plot: Raw Het
# -----------------------------
ax_raw.scatter(positions, female_het, color=female_color, label="Females (raw)", alpha=0.5, s=8)
ax_raw.scatter(positions, male_het, color=male_color, label="Males (raw)", alpha=0.5, s=8)
ax_raw.axvspan(args.region_start, args.region_end, color=region_color, alpha=0.1)
ax_raw.set_xlim(xlim_raw)
ax_raw.set_ylabel("Raw Het")
ax_raw.set_title("Raw Heterozygosity")
ax_raw.tick_params(labelbottom=False)

# -----------------------------
# Plot: Gene Region
# -----------------------------
for _, row in region_sorted.iterrows():
    y = row["track"] * 0.5
    start, end, strand = row["start"], row["end"], row["strand"]
    label = row["gene_label"][:20]
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
    ax_gene.add_patch(arrow)
    ax_gene.text((start + end) / 2, y + 0.2, label, ha="center", va="bottom", fontsize=7)

ax_gene.set_ylim(-0.5, max(region_sorted["track"]) * 0.6 + 1)
ax_gene.set_yticks([])
ax_gene.set_ylabel("Genes")
ax_gene.set_xlim(xlim_full)
ax_gene.tick_params(labelbottom=False)

# -----------------------------
# Plot: Smoothed Het
# -----------------------------
ax_smooth.plot(positions, female_het_smooth, label="Females (smoothed)", color=female_color)
ax_smooth.plot(positions, male_het_smooth, label="Males (smoothed)", color=male_color)
ax_smooth.plot(positions, diff_het, label="Females - Males", color=diff_color, linestyle="--")
ax_smooth.axvspan(args.region_start, args.region_end, color=region_color, alpha=0.2)
ax_smooth.set_xlabel("Genomic Position (bp)")
ax_smooth.set_ylabel("Smoothed Het")
ax_smooth.set_title("Smoothed Heterozygosity")
ax_smooth.set_xlim(xlim_full)

# -----------------------------
# Legend
# -----------------------------
legend_elements = [
    mpatches.Patch(color=female_color, label="Females"),
    mpatches.Patch(color=male_color, label="Males"),
    mpatches.Patch(color=diff_color, label="Females - Males"),
    mpatches.Patch(color=arrow_color, label="Gene Arrow"),
    mpatches.Patch(color=region_color, label="Region of Interest")
]
ax_legend.legend(handles=legend_elements, loc="center", frameon=False, fontsize=10)

# -----------------------------
# Save
# -----------------------------
plt.tight_layout()
plt.savefig(args.out_png)
plt.savefig(args.out_pdf)
plt.savefig(args.out_svg)
