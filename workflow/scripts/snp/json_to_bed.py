import json
import argparse

def translate_chrom_name(chrom):
    mapping = {
        "chr4": "OV696689.1"
    }
    return mapping.get(chrom, chrom)

def main(json_file, bed_file):
    with open(json_file) as f:
        snps = json.load(f)

    with open(bed_file, "w") as out:
        for snp in snps:
            chrom = translate_chrom_name(snp["chromosome"])
            start = int(snp["position"]) - 1  # BED is 0-based
            end = int(snp["position"])
            pval = snp.get("p_value", ".")
            out.write(f"{chrom}\t{start}\t{end}\t{pval}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--json", required=True)
    parser.add_argument("--bed", required=True)
    args = parser.parse_args()
    main(args.json, args.bed)
