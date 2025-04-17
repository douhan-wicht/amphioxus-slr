import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True, help="Locus-tag based gene list")
parser.add_argument("--map", required=True, help="locus_tag → GeneID mapping file")
parser.add_argument("--output", required=True)
args = parser.parse_args()

# Read locus_tag → GeneID mapping
mapping = pd.read_csv(args.map, sep="\t", header=None, names=["locus_tag", "geneid"])
map_dict = dict(zip(mapping.locus_tag, mapping.geneid.astype(int)))

# Read gene list and convert
with open(args.input) as f:
    genes = [line.strip().replace("gene-", "") for line in f if line.strip()]

converted = [map_dict[gene] for gene in genes if gene in map_dict]

with open(args.output, "w") as out:
    for geneid in converted:
        out.write(str(geneid) + "\n")
