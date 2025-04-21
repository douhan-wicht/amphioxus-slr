import argparse
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

def parse_args():
    parser = argparse.ArgumentParser(description="Classify SNP effects manually.")
    parser.add_argument("--snps", required=True, help="Path to SNPs TSV file with CHROM, POS, REF, ALT columns")
    parser.add_argument("--gff", required=True, help="Path to GFF file with CDS annotations")
    parser.add_argument("--ref", required=True, help="Path to reference genome in FASTA format")
    parser.add_argument("--output", required=True, help="Output TSV file with SNP effects")
    return parser.parse_args()

# Map CHROM values like "chr4" to actual FASTA IDs
chromosome_map = {
    "chr4": "OV696689.1"
}

def classify_snp(chrom, pos, ref_base, alt_base, ref_genome, cds):
    pos = int(pos)
    fasta_id = chromosome_map.get(chrom, chrom)
    if fasta_id not in ref_genome:
      raise KeyError(f"Chromosome {fasta_id} not found in reference genome.")
    ref_seq = ref_genome[fasta_id].seq


    overlapping = cds[(cds["seqid"] == chrom) & (cds["start"] <= pos) & (cds["end"] >= pos)]
    if overlapping.empty:
        return "intergenic"

    row = overlapping.iloc[0]
    strand = row["strand"]
    phase = int(row["phase"]) if row["phase"] != '.' else 0
    cds_start = row["start"]

    rel_pos = pos - cds_start
    codon_start = pos - ((rel_pos + phase) % 3)
    codon_seq = ref_seq[codon_start:codon_start + 3]

    if len(codon_seq) != 3:
        return "partial"

    if strand == "-":
        codon_seq = codon_seq.reverse_complement()
        snp_index = 2 - ((pos - codon_start) % 3)
    else:
        snp_index = (pos - codon_start) % 3

    ref_codon = str(codon_seq)
    alt_codon = list(ref_codon)
    alt_codon[snp_index] = alt_base
    alt_codon = ''.join(alt_codon)

    aa_ref = str(Seq(ref_codon).translate())
    aa_alt = str(Seq(alt_codon).translate())

    if aa_ref == aa_alt:
        return "silent"
    elif aa_alt == "*":
        return "nonsense"
    else:
        return "missense"

def main():
    args = parse_args()

    snps = pd.read_csv(args.snps, sep="\t", usecols=["CHROM", "POS", "REF", "ALT"])
    gff = pd.read_csv(args.gff, sep="\t", comment='#', header=None,
                      names=["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"])
    cds = gff[gff["type"] == "CDS"]

    ref_genome = SeqIO.to_dict(SeqIO.parse(args.ref, "fasta"))

    snps["effect"] = snps.apply(lambda row: classify_snp(
        row["CHROM"], row["POS"], row["REF"], row["ALT"], ref_genome, cds), axis=1)

    snps.to_csv(args.output, sep="\t", index=False)

if __name__ == "__main__":
    main()