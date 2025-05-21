#!/usr/bin/env python

import argparse
import pandas as pd
import glob
import os

parser = argparse.ArgumentParser()
parser.add_argument("--input", nargs="+", required=True, help="List of input count files")
parser.add_argument("--output", required=True, help="Path to output summary table")
args = parser.parse_args()

data = []
for filepath in args.input:
    filename = os.path.basename(filepath)
    chrom = filename.replace("pass_snps_", "").replace(".txt", "")
    
    with open(filepath) as f:
        count = int(f.read().strip())
    data.append((chrom, count))

df = pd.DataFrame(data, columns=["Chromosome", "PASS_SNP_Count"]).sort_values(by="Chromosome")
df.to_csv(args.output, sep="\t", index=False)
