import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import matplotlib.patches as mpatches

# -----------------------------
# STYLE SETUP
# -----------------------------
sns.set(style="whitegrid", context="talk", palette="colorblind")
colors = sns.color_palette("colorblind")
female_color = colors[0]   # blue-ish
male_color = colors[2]     # green-ish
diff_color = colors[4]     # purple-ish
region_color = colors[3]   # red-ish

# -----------------------------
# ARGPARSE
# -----------------------------
parser = argparse.ArgumentParser(description="Combine smoothed heterozygosity plot with gene annotation track")

parser.add_argument("--vcf_tab", required=True, help="VCF tabular file with genotypes")
parser.add_argument("--gff", required=True, help="Path to GFF3 annotation file")
parser.add_argument("--seqid", required=True, help="Chromosome/scaffold name (e.g., 'chr4')")
parser.add_argument("--region_start", type=int, required=True, help="Start coordinate of region")
parser.add_argument("--region_end", type=int, required=True, help="End coordinate of region")
parser.add_argument("--out_png", required=True, help="Output PNG path")
parser.add_argument("--out_pdf", required=True, help="Output PDF path")

args = parser.parse_args()

# -----------------------------
# LOAD DATA
# -----------------------------
df = pd.read_csv(args.vcf_tab, sep="\t")

col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
gff = pd.read_csv(args.gff, sep="\t", comment="#", names=col_names)

# -----------------------------
# HELPER FUNCTIONS
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

def extract_gene_name(attr):
    for field in attr.split(";"):
        if "Name=" in field:
            return field.split("Name=")[-1]
        elif "gene=" in field:
            return field.split("gene=")[-1]
        elif "ID=" in field:
            return field.split("ID=")[-1]
    return "unknown"

# -----------------------------
# HETEROZYGOSITY
# -----------------------------
gt_columns = [col for col in df.columns if col.endswith(".GT")]
sex_map = {col: infer_sex(col) for col in gt_columns}
female_inds = [col for col in gt_columns if sex_map[col] == "F"]
male_inds = [col for col in gt_columns if sex_map[col] == "M"]

positions = df["POS"]
female_het = df[female_inds].applymap(compute_heterozygosity).mean(axis=1)
male_het = df[male_inds].applymap(compute_heterozygosity).mean(axis=1)

# Smooth
window_size = 50
female_het_smooth = female_het.rolling(window=window_size, min_periods=1).mean()
male_het_smooth = male_het.rolling(window=window_size, min_periods=1).mean()
diff_het = female_het_smooth - male_het_smooth

# -----------------------------
# GENE ANNOTATIONS
# -----------------------------
valid_types = ["gene", "mRNA", "transcript"]
region = gff[
    (gff["seqid"] == args.seqid) &
    (gff["type"].isin(valid_types)) &
    (gff["start"] <= args.region_end) &
    (gff["end"] >= args.region_start)
].copy()

region["gene_name"] = region["attributes"].apply(extract_gene_name)
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
# PLOTTING
# -----------------------------
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 9), height_ratios=[3, 1], sharex=True)

# --- Heterozygosity Plot ---
ax1.plot(positions, female_het_smooth, label="Females (smoothed)", color=female_color)
ax1.plot(positions, male_het_smooth, label="Males (smoothed)", color=male_color)
ax1.plot(positions, diff_het, label="Females - Males", color=diff_color, linestyle="--")
ax1.axvspan(args.region_start, args.region_end, color=region_color, alpha=0.2, label="Region of Interest")

ax1.set_ylabel("Mean Heterozygosity")
ax1.set_title("Smoothed Heterozygosity by Sex Across Genomic Region")
ax1.legend()
ax1.grid(True)

# --- Gene Annotation Track ---
ax2.axvspan(args.region_start, args.region_end, color=region_color, alpha=0.2)

for _, row in region_sorted.iterrows():
    y = row["track"] * 0.5
    strand = row["strand"]
    color = female_color if strand == "+" else male_color  # arbitrary use of distinct palette colors
    name = row["gene_name"][:20]
    start = row["start"]
    end = row["end"]

    arrow = mpatches.FancyArrow(
        start if strand == "+" else end,
        y,
        end - start if strand == "+" else start - end,
        0,
        width=0.1,
        length_includes_head=True,
        head_width=0.2,
        head_length=500,
        color=color,
        alpha=0.9
    )
    ax2.add_patch(arrow)
    ax2.text((start + end) / 2, y + 0.2, name, ha="center", va="bottom", fontsize=7, rotation=0)

ax2.set_ylim(-0.5, max(region_sorted["track"]) * 0.6 + 1)
ax2.set_yticks([])
ax2.set_xlabel("Base Pair Position (bp)")
ax2.set_title(f"Genes in {args.seqid}:{args.region_start}-{args.region_end}")

# Shared X range
x_min = args.region_start - 5000
x_max = args.region_end + 100000
ax1.set_xlim(x_min, x_max)
ax2.set_xlim(x_min, x_max)

plt.tight_layout()
plt.savefig(args.out_png)
plt.savefig(args.out_pdf)
