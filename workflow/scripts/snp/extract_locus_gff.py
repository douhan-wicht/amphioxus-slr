import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--loci", nargs="+", required=True)
    args = parser.parse_args()

    loci_set = set(args.loci)

    kept_lines = []
    with open(args.input) as fin:
        for line in fin:
            if line.startswith("#"):
                kept_lines.append(line)
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue
            attr = fields[8]
            if any(locus in attr for locus in loci_set):
                kept_lines.append(line)

    with open(args.output, "w") as fout:
        fout.writelines(kept_lines)

    print(f"Wrote {len(kept_lines)} lines for loci: {', '.join(args.loci)}")

if __name__ == "__main__":
    main()
