import argparse
import pandas as pd
import json

def parse_args():
    parser = argparse.ArgumentParser(description="Extract top SNPs from master list.")
    parser.add_argument("--json", required=True, help="Path to JSON with top SNPs")
    parser.add_argument("--snps", required=True, help="Path to full SNP TSV (must have CHROM, POS, REF, ALT)")
    parser.add_argument("--output", required=True, help="Path to output TSV of filtered top SNPs")
    return parser.parse_args()

def main():
    args = parse_args()

    # Load top SNPs positions from JSON
    with open(args.json) as f:
        top_snps = json.load(f)
    top_positions = {(snp["chromosome"], snp["position"]) for snp in top_snps}

    # Load full SNP table
    df = pd.read_csv(args.snps, sep="\t")

    # Filter
    filtered = df[df.apply(lambda row: (row["CHROM"], row["POS"]) in top_positions, axis=1)]

    # Output only columns required for annotation
    columns = ["CHROM", "POS", "REF", "ALT"]
    filtered[columns].to_csv(args.output, sep="\t", index=False)

if __name__ == "__main__":
    main()