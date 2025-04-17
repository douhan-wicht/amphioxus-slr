import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--mapping", required=True, help="GFF extracted locus_tag ↔ gene ID")
parser.add_argument("--gene_info", required=True, help="NCBI filtered gene_info file")
parser.add_argument("--output", required=True, help="Output mapping file")
args = parser.parse_args()

# GFF geneid-to-locus tag mapping
df_gff = pd.read_csv(args.mapping, sep="\t", header=None, names=["GeneID", "locus_tag"])

# NCBI gene_info file (filtered)
columns = [
    "tax_id", "GeneID", "Symbol", "LocusTag", "Synonyms", "dbXrefs",
    "chromosome", "map_location", "description", "type_of_gene",
    "Symbol_from_nomenclature_authority", "Full_name_from_nomenclature_authority",
    "Nomenclature_status", "Other_designations", "Modification_date", "Feature_type"
]
df_info = pd.read_csv(args.gene_info, sep="\t", comment="#", names=columns)

print("GFF head:")
print(df_gff.head())
print("gene_info head:")
print(df_info.head())


# Merge using locus tag
df_merge = pd.merge(df_gff, df_info, left_on="locus_tag", right_on="LocusTag")

# Output: final map from locus_tag to GeneID
df_merge[["locus_tag", "GeneID"]].drop_duplicates().to_csv(args.output, sep="\t", index=False, header=False)
