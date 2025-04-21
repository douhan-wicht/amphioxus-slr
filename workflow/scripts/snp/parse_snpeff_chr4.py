import argparse
import vcf  # pyvcf
import csv
from collections import defaultdict

# ⬇️ Parse CLI args passed by Snakemake shell
parser = argparse.ArgumentParser()
parser.add_argument('--input', required=True, help="SnpEff-annotated VCF file")
parser.add_argument('--output', required=True, help="Output TSV file")
args = parser.parse_args()

vcf_file = args.input
out_file = args.output

LOCUS_TO_GENE = {
    "BLAG_LOCUS17194": "HAO1",
    "BLAG_LOCUS17195": "FLT1"
}

# To hold filtered variants
summary_data = []

# Read the annotated VCF
vcf_reader = vcf.Reader(filename=vcf_file)

for record in vcf_reader:
    if not record.INFO.get("ANN"):
        continue

    annotations = record.INFO["ANN"]
    for ann in annotations:
        fields = ann.split('|')
        effect = fields[1]
        impact = fields[2]
        gene_id = fields[3]  # This is like BLAG_LOCUS1234
        feature_id = fields[4]
        feature_type = fields[5]
        exon_rank = fields[8]
        cdna_change = fields[9]
        aa_change = fields[10]
        protein_pos = fields[13]

        # ⬇️ Check if this gene is one of interest
        if gene_id not in LOCUS_TO_GENE:
            continue

        gene_name = LOCUS_TO_GENE[gene_id]

        summary_data.append({
            "CHROM": record.CHROM,
            "POS": record.POS,
            "REF": record.REF,
            "ALT": ','.join(str(alt) for alt in record.ALT),
            "GENE": gene_name,
            "IMPACT": impact,
            "EFFECT": effect,
            "EXON_RANK": exon_rank,
            "PROTEIN_CHANGE": aa_change,
            "CDNA_CHANGE": cdna_change,
            "PROTEIN_POS": protein_pos
        })

# ⬇️ Write output only if we found matching annotations
if summary_data:
    with open(out_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=summary_data[0].keys(), delimiter='\t')
        writer.writeheader()
        writer.writerows(summary_data)
else:
    print("No annotations found for target genes.")
