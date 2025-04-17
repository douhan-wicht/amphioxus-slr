# scripts/go/convert_gene_list_to_ids.py
import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--gene_list", required=True)
parser.add_argument("--mapping", required=True)
parser.add_argument("--output", required=True)
args = parser.parse_args()

# Load list of gene IDs (e.g. gene-BLAG_LOCUS####)
with open(args.gene_list) as f:
    gene_list = set(line.strip() for line in f)

# Load mapping (gene_id → GeneID)
df_map = pd.read_csv(args.mapping, sep="\t", header=None, names=["gene_id", "GeneID"])

# Filter and convert
df_filtered = df_map[df_map["gene_id"].isin(gene_list)]
df_filtered["GeneID"].astype(int).to_csv(args.output, index=False, header=False)
