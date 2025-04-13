import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# -----------------------------
# ARGPARSE
# -----------------------------
parser = argparse.ArgumentParser(description="Plot smoothed heterozygosity by sex")

parser.add_argument("--input", required=True, help="Input .tab file with genotypes")
parser.add_argument("--out_png", required=True, help="Output PNG path")
parser.add_argument("--out_pdf", required=True, help="Output PDF path")
parser.add_argument("--region_start", type=int, default=6142346, help="Start of region of interest")
parser.add_argument("--region_end", type=int, default=6164195, help="End of region of interest")

args = parser.parse_args()

# -----------------------------
# LOAD DATA
# -----------------------------
df = pd.read_csv(args.input, sep="\t")

# Function to compute heterozygosity
def compute_heterozygosity(gt):
    if not isinstance(gt, str) or gt in ("./.", ".|."):
        return np.nan
    alleles = gt.replace("|", "/").split("/")
    if len(alleles) != 2:
        return np.nan
    return 1 if alleles[0] != alleles[1] else 0

# Detect genotype columns
gt_columns = [col for col in df.columns if col.endswith(".GT")]

# Infer sex from sample names
def infer_sex(name):
    name = name.upper()
    if name.startswith("F") or name.startswith("RF"):
        return "F"
    elif name.startswith("M") or name.startswith("RM"):
        return "M"
    else:
        return "U"

sex_map = {col: infer_sex(col) for col in gt_columns}
female_inds = [col for col in gt_columns if sex_map[col] == "F"]
male_inds = [col for col in gt_columns if sex_map[col] == "M"]

# Compute heterozygosity
positions = df["POS"]
female_het = df[female_inds].applymap(compute_heterozygosity).mean(axis=1)
male_het = df[male_inds].applymap(compute_heterozygosity).mean(axis=1)

# Smooth
window_size = 50
female_het_smooth = female_het.rolling(window=window_size, min_periods=1).mean()
male_het_smooth = male_het.rolling(window=window_size, min_periods=1).mean()
diff_het = female_het_smooth - male_het_smooth

# Plot
plt.figure(figsize=(14, 7))
plt.plot(positions, female_het_smooth, label="Females (smoothed)", color="blue")
plt.plot(positions, male_het_smooth, label="Males (smoothed)", color="green")
plt.plot(positions, diff_het, label="Females - Males", color="purple", linestyle="--")

plt.axvspan(args.region_start, args.region_end, color="red", alpha=0.2, label="Region of Interest")

plt.xlabel("Base Pair Position")
plt.ylabel("Mean Heterozygosity")
plt.title("Smoothed Heterozygosity by Sex Across Genomic Region")
plt.grid(True)
plt.legend()
plt.xlim(args.region_start - 5000, args.region_end + 5000)
plt.tight_layout()
plt.savefig(args.out_png)
plt.savefig(args.out_pdf)
