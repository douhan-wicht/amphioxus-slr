import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# -----------------------------
# COLORBLIND-FRIENDLY SETTINGS
# -----------------------------
sns.set(style="whitegrid", context="talk", palette="colorblind")

# -----------------------------
# ARGPARSE
# -----------------------------
parser = argparse.ArgumentParser(description="Plot raw heterozygosity by sex (no smoothing)")
parser.add_argument("--input", required=True, help="Input .tab file with genotypes")
parser.add_argument("--out_png", required=True, help="Output PNG path")
parser.add_argument("--out_pdf", required=True, help="Output PDF path")
parser.add_argument("--out_svg", required=True, help="Output SVG path")
parser.add_argument("--region_start", type=int, default=6142346, help="Start of region of interest")
parser.add_argument("--region_end", type=int, default=6164195, help="End of region of interest")
args = parser.parse_args()

# -----------------------------
# LOAD DATA
# -----------------------------
df = pd.read_csv(args.input, sep="\t")

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

# -----------------------------
# PROCESSING
# -----------------------------
gt_columns = [col for col in df.columns if col.endswith(".GT")]
sex_map = {col: infer_sex(col) for col in gt_columns}
female_inds = [col for col in gt_columns if sex_map[col] == "F"]
male_inds = [col for col in gt_columns if sex_map[col] == "M"]
positions = df["POS"]

female_het = df[female_inds].applymap(compute_heterozygosity).mean(axis=1)
male_het = df[male_inds].applymap(compute_heterozygosity).mean(axis=1)

# -----------------------------
# PLOTTING
# -----------------------------
cb_palette = sns.color_palette("colorblind")
female_color = cb_palette[0]
male_color = cb_palette[3]

plt.figure(figsize=(14, 6))
plt.scatter(positions, female_het, label="Females", color=female_color, alpha=0.6, s=8)
plt.scatter(positions, male_het, label="Males", color=male_color, alpha=0.6, s=8)
plt.axvspan(args.region_start, args.region_end, color=cb_palette[6], alpha=0.1, label="Region of Interest")

plt.xlabel("Base Pair Position")
plt.ylabel("Mean Heterozygosity")
plt.title("Raw Mean Heterozygosity by Sex")
plt.grid(True)
plt.xlim(args.region_start - 5000, args.region_end + 5000)
plt.tight_layout()
plt.savefig(args.out_png)
plt.savefig(args.out_pdf)
plt.savefig(args.out_svg)
