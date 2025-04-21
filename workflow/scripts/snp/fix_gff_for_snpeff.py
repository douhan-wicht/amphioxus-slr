#!/usr/bin/env python3

import argparse
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description="Fix GFF for SnpEff.")
    parser.add_argument("--input", dest="input_gff", type=Path, required=True, help="Input GFF file")
    parser.add_argument("--output", dest="output_gff", type=Path, required=True, help="Output GFF file")
    args = parser.parse_args()

    gene_to_rna = {}

    with args.input_gff.open() as infile, args.output_gff.open("w") as outfile:
        for line in infile:
            if line.startswith("#") or not line.strip():
                outfile.write(line)
                continue

            fields = line.strip().split("\t")
            if len(fields) != 9:
                outfile.write(line)
                continue

            attr_field = fields[8]
            attrs = dict(x.split("=", 1) for x in attr_field.split(";") if "=" in x)

            if fields[2] == "mRNA":
                rna_id = attrs.get("ID")
                gene_id = attrs.get("Parent")
                gene_to_rna[gene_id] = rna_id

                if "Name" not in attrs and "gene" in attrs:
                    attrs["Name"] = attrs["gene"]

            elif fields[2] in {"CDS", "exon", "five_prime_UTR", "three_prime_UTR"}:
                parent = attrs.get("Parent")
                if parent in gene_to_rna:
                    attrs["Parent"] = gene_to_rna[parent]

            fields[8] = ";".join(f"{k}={v}" for k, v in attrs.items())
            outfile.write("\t".join(fields) + "\n")

if __name__ == "__main__":
    main()
