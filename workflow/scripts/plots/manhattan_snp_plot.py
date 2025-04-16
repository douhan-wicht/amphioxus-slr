import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import fisher_exact

# -----------------------------
# STYLE: Colorblind Palette
# -----------------------------
sns.set(style="whitegrid", context="talk", palette="colorblind")
colors = sns.color_palette("colorblind", n_colors=20)

# -----------------------------
# ARGPARSE
# -----------------------------
parser = argparse.ArgumentParser(description="Genome-wide Manhattan plot of SNP sex-association p-values")
parser.add_argument("--input", required=True, help="Input .tab file with genotypes")
parser.add_argument("--out_png", required=True, help="Output PNG path")
parser.add_argument("--out_pdf", required=True, help="Output PDF path")
args = parser.parse_args()

# -----------------------------
# LOAD DATA
# -----------------------------
df = pd.read_csv(args.input, sep="\t")
gt_columns = [col for col in df.columns if col.endswith(".GT")]

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

def parse_alleles(geno):
    if not isinstance(geno, str) or geno in ["./.", ".|.", ".", ""]:
        return None
    return geno.replace("|", "/").split("/")

def has_alt(geno, alt):
    alleles = parse_alleles(geno)
    return alleles is not None and alt in alleles

# -----------------------------
# COMPUTE P-VALUES
# -----------------------------
results = []

for idx, row in df.iterrows():
    alt = row["ALT"]
    male_alt = male_ref = female_alt = female_ref = 0

    for col in female_inds:
        if has_alt(row[col], alt):
            female_alt += 1
        elif parse_alleles(row[col]) is not None:
            female_ref += 1

    for col in male_inds:
        if has_alt(row[col], alt):
            male_alt += 1
        elif parse_alleles(row[col]) is not None:
            male_ref += 1

    table = [[female_ref, female_alt], [male_ref, male_alt]]
    if all(sum(row) > 0 for row in table):
        _, pval = fisher_exact(table)
    else:
        pval = np.nan

    results.append({
        "chr": row["CHROM"],
        "pos": row["POS"],
        "pval": pval
    })

res_df = pd.DataFrame(results).dropna()
res_df["-log10p"] = -np.log10(res_df["pval"])

# -----------------------------
# COMPUTE CUMULATIVE POSITIONS
# -----------------------------
res_df["chr_num"] = res_df["chr"].str.replace("chr", "").astype(int)
res_df = res_df.sort_values(["chr_num", "pos"])

chromosomes = res_df["chr"].unique()
chrom_offsets = {}
offset = 0
cumulative_pos = []

for chrom in sorted(chromosomes, key=lambda x: int(x.replace("chr", ""))):
    chr_df = res_df[res_df["chr"] == chrom]
    chrom_offsets[chrom] = offset
    cumulative_pos.extend(chr_df["pos"] + offset)
    offset += chr_df["pos"].max() + 1e6  # padding

res_df["cumulative_pos"] = cumulative_pos

# -----------------------------
# PLOTTING
# -----------------------------
fig, ax = plt.subplots(figsize=(16, 7))

for i, chrom in enumerate(sorted(res_df["chr"].unique(), key=lambda x: int(x.replace("chr", "")))):
    chrom_data = res_df[res_df["chr"] == chrom]
    ax.scatter(
        chrom_data["cumulative_pos"],
        chrom_data["-log10p"],
        color=colors[i % len(colors)],
        s=3,
        label=chrom if i % 2 == 0 else "",  # only label every second chr
        alpha=0.7
    )

# Threshold line
ax.axhline(-np.log10(5e-8), color="darkorange", linestyle="--", label="Genome-wide threshold")

# Ticks and labels
xticks = []
xlabels = []

for chrom in sorted(res_df["chr"].unique(), key=lambda x: int(x.replace("chr", ""))):
    chr_df = res_df[res_df["chr"] == chrom]
    center = chr_df["cumulative_pos"].median()
    xticks.append(center)
    xlabels.append(chrom)

ax.set_xticks(xticks)
ax.set_xticklabels(xlabels, rotation=45)
ax.set_xlabel("Genomic Position (by Chromosome)")
ax.set_ylabel("-log10(p-value)")
ax.set_title("Genome-wide Sex-Associated SNPs", weight="bold")
ax.legend()
ax.grid(True, axis="y")

plt.tight_layout()
plt.savefig(args.out_png)
plt.savefig(args.out_pdf)
