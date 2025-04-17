# scripts/go/gff_extract_locus_tags.py
import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--gff", required=True)
parser.add_argument("--output", required=True)
args = parser.parse_args()

rows = []
with open(args.gff) as f:
    for line in f:
        if line.startswith("#"):
            continue
        fields = line.strip().split("\t")
        if len(fields) < 9 or fields[2] != "gene":
            continue
        attr_field = fields[8]
        attrs = dict(part.split("=", 1) for part in attr_field.split(";") if "=" in part)
        if "ID" in attrs and "locus_tag" in attrs:
            rows.append((attrs["ID"], attrs["locus_tag"]))

df = pd.DataFrame(rows, columns=["gene_id", "locus_tag"])
df.to_csv(args.output, sep="\t", index=False, header=False)
