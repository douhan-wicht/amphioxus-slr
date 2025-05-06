import pandas as pd
import argparse

def classify_genotype(gt):
    if gt in ["./.", ".", ""]:
        return "missing"
    alleles = gt.replace("|", "/").split("/")
    return "heterozygous" if alleles[0] != alleles[1] else "homozygous"

def analyze(file, start, end, pos_list, output_prefix):
    df = pd.read_csv(file, sep="\t")

    # Optional filtering by region
    if start is not None and end is not None:
        df = df[(df["POS"] >= start) & (df["POS"] <= end)]

    # Optional filtering by SNP list
    if pos_list is not None:
        snp_positions = pd.read_csv(pos_list, header=None)[0].tolist()
        df = df[df["POS"].isin(snp_positions)]

    # Only use known samples
    female_cols = [col for col in df.columns if (col.startswith("F") or col.startswith("RF")) and not col.startswith("RU")]
    male_cols = [col for col in df.columns if (col.startswith("M") or col.startswith("RM")) and not col.startswith("RU")]

    female_class = df[female_cols].applymap(classify_genotype)
    male_class = df[male_cols].applymap(classify_genotype)

    female_summary = female_class.apply(pd.Series.value_counts, axis=1).fillna(0)
    male_summary = male_class.apply(pd.Series.value_counts, axis=1).fillna(0)

    result = df[["CHROM", "POS"]].copy()
    result["female_heterozygous"] = female_summary.get("heterozygous", 0)
    result["female_homozygous"] = female_summary.get("homozygous", 0)
    result["male_heterozygous"] = male_summary.get("heterozygous", 0)
    result["male_homozygous"] = male_summary.get("homozygous", 0)

    result.to_csv(f"{output_prefix}.tsv", sep="\t", index=False)

    # Summary stats
    snp_count = len(result)
    avg_female_het = result["female_heterozygous"].mean()
    avg_male_hom = result["male_homozygous"].mean()
    prop_female_het_10 = (result["female_heterozygous"] >= 10).mean() * 100
    prop_male_hom_15 = (result["male_homozygous"] >= 15).mean() * 100

    summary = f"""
    Filtered by:
        Region: {start} - {end} (if provided)
        SNP List: {pos_list if pos_list else 'None'}

    Total SNPs after filtering: {snp_count}

    Female (heterozygous dominant):
        - Avg. heterozygous per SNP: {avg_female_het:.2f}
        - % SNPs with >=10 heterozygous females: {prop_female_het_10:.1f}%

    Male (homozygous dominant):
        - Avg. homozygous per SNP: {avg_male_hom:.2f}
        - % SNPs with >=15 homozygous males: {prop_male_hom_15:.1f}%
    """
    with open(f"{output_prefix}_summary.txt", "w") as f:
        f.write(summary.strip())

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--start", type=int, default=None)
    parser.add_argument("--end", type=int, default=None)
    parser.add_argument("--pos-list", default=None)
    parser.add_argument("--output-prefix", required=True)
    args = parser.parse_args()

    analyze(args.input, args.start, args.end, args.pos_list, args.output_prefix)
