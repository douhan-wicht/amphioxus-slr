import pandas as pd
from goatools.obo_parser import GODag
from goatools.goea.go_enrichment_ns import GOEnrichmentStudy
from goatools.associations import read_ncbi_gene2go
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--gene2go", required=True)
parser.add_argument("--target", required=True)
parser.add_argument("--background", required=True)
parser.add_argument("--obo", required=True)
parser.add_argument("--out", required=True)
args = parser.parse_args()

# Read inputs
gene2go = read_ncbi_gene2go(args.gene2go, taxids=[7740])
go_dag = GODag(args.obo)

with open(args.background) as f:
    background_genes = [int(line.strip()) for line in f if line.strip().isdigit()]

with open(args.target) as f:
    target_genes = [int(line.strip()) for line in f if line.strip().isdigit()]

# Enrichment
goeaobj = GOEnrichmentStudy(
    background_genes, gene2go, go_dag,
    propagate_counts=True,
    alpha=0.05,
    methods=['fdr_bh']
)

results = goeaobj.run_study(target_genes)

# Save results
df = pd.DataFrame([
    {
        "GO": r.GO,
        "NS": r.NS,
        "name": r.name,
        "p_fdr_bh": r.p_fdr_bh,
        "ratio_in_study": f"{r.study_count}/{r.study_n}",
        "ratio_in_pop": f"{r.pop_count}/{r.pop_n}",
        "enrichment": r.enrichment,
    }
    for r in results if r.p_fdr_bh < 0.05
])
df.sort_values("p_fdr_bh", inplace=True)
df.to_csv(args.out, sep="\t", index=False)
