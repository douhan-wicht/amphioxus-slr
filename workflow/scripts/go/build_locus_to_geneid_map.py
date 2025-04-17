import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--gff", required=True, help="GFF-derived gene ID to locus_tag mapping")
parser.add_argument("--gene_info", required=True, help="Filtered gene_info.tsv")
parser.add_argument("--output", required=True, help="Output TSV mapping locus_tag to GeneID")
args = parser.parse_args()

# Step 1: Read GFF mapping file (auto-detect delimiters, fallback to whitespace)
with open(args.gff) as f:
    lines = [line.strip().split()[:2] for line in f if line.strip() and not line.startswith("#")]

df_gff = pd.DataFrame(lines, columns=["GeneID", "locus_tag"])

# Step 2: Read gene_info file with header
df_info = pd.read_csv(args.gene_info, sep="\t")

# Step 3: Merge on locus_tag
df_merge = pd.merge(df_gff, df_info, on="locus_tag", how="inner")

# Step 4: Output
df_merge[["locus_tag", "GeneID_y"]].drop_duplicates().to_csv(args.output, sep="\t", index=False, header=False)
