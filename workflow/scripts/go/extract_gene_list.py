import argparse
import pandas as pd

# -----------------------------
# ARGPARSE
# -----------------------------
parser = argparse.ArgumentParser(description="Extract gene list from GFF3 in a specific region")
parser.add_argument("--gff", required=True, help="Path to GFF3 file")
parser.add_argument("--seqid", required=True, help="Chromosome/scaffold ID (e.g., OV696689.1)")
parser.add_argument("--start", type=int, required=True, help="Region start coordinate")
parser.add_argument("--end", type=int, required=True, help="Region end coordinate")
parser.add_argument("--output", required=True, help="Output file with gene list")

args = parser.parse_args()

# -----------------------------
# Load GFF
# -----------------------------
col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
gff = pd.read_csv(args.gff, sep="\t", comment="#", names=col_names)

# -----------------------------
# Filter genes in region
# -----------------------------
region_genes = gff[
    (gff["seqid"] == args.seqid) &
    (gff["type"] == "gene") &
    (gff["start"] <= args.end) &
    (gff["end"] >= args.start)
].copy()

# -----------------------------
# Extract gene names
# -----------------------------
def extract_gene_name(attr):
    for field in attr.split(";"):
        if "Name=" in field:
            return field.split("Name=")[-1]
        elif "gene=" in field:
            return field.split("gene=")[-1]
        elif "ID=" in field:
            return field.split("ID=")[-1]
    return "unknown"

region_genes["gene_name"] = region_genes["attributes"].apply(extract_gene_name)

# -----------------------------
# Save output
# -----------------------------
region_genes["gene_name"].dropna().drop_duplicates().to_csv(args.output, index=False, header=False)
print(f"✅ Extracted {len(region_genes)} genes to {args.output}")
