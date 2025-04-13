import argparse
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.patches as mpatches

# -----------------------------
# ARGPARSE
# -----------------------------
parser = argparse.ArgumentParser(description="Plot genes in a region from GFF3")

parser.add_argument("--gff", required=True, help="Path to GFF3 annotation file")
parser.add_argument("--seqid", required=True, help="Chromosome/scaffold name (e.g., 'chr4')")
parser.add_argument("--start", type=int, required=True, help="Start coordinate of region")
parser.add_argument("--end", type=int, required=True, help="End coordinate of region")
parser.add_argument("--out_png", required=True, help="Output PNG path")
parser.add_argument("--out_pdf", required=True, help="Output PDF path")

args = parser.parse_args()

# -----------------------------
# LOAD AND INSPECT GFF
# -----------------------------
col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
gff = pd.read_csv(args.gff, sep="\t", comment="#", names=col_names)

# Debug output
print("\n🧬 DEBUG INFO")
print("Unique seqids in GFF:", gff["seqid"].unique())
print(f"Feature count on {args.seqid}:", (gff["seqid"] == args.seqid).sum())
print(f"Feature types on {args.seqid}:\n", gff[gff["seqid"] == args.seqid]["type"].value_counts())
print(f"First few entries on {args.seqid}:\n", gff[gff["seqid"] == args.seqid].head())

# -----------------------------
# FILTER REGION & FEATURES
# -----------------------------
valid_types = ["gene", "mRNA", "transcript"]
region = gff[
    (gff["seqid"] == args.seqid) &
    (gff["type"].isin(valid_types)) &
    (gff["start"] <= args.end) &
    (gff["end"] >= args.start)
].copy()

print(f"\n🧬 Number of features in region {args.seqid}:{args.start}-{args.end}:", len(region))

# -----------------------------
# PARSE GENE NAMES
# -----------------------------
def extract_gene_name(attr):
    for field in attr.split(";"):
        if "Name=" in field:
            return field.split("Name=")[-1]
        elif "gene=" in field:
            return field.split("gene=")[-1]
        elif "ID=" in field:
            return field.split("ID=")[-1]
    return "unknown"

region["gene_name"] = region["attributes"].apply(extract_gene_name)

# -----------------------------
# PLOTTING (improved)
# -----------------------------
fig, ax = plt.subplots(figsize=(12, 2.5))

# Sort genes by start and assign "tracks" (stacking rows)
region_sorted = region.sort_values("start").reset_index(drop=True)
region_sorted["track"] = 0
track_ends = []

for i, row in region_sorted.iterrows():
    placed = False
    for track_idx, end_pos in enumerate(track_ends):
        if row["start"] > end_pos + 1000:  # gap between genes
            region_sorted.at[i, "track"] = track_idx
            track_ends[track_idx] = row["end"]
            placed = True
            break
    if not placed:
        region_sorted.at[i, "track"] = len(track_ends)
        track_ends.append(row["end"])

# Plot baseline
ax.hlines(y=0, xmin=args.start, xmax=args.end, color="black", linewidth=2)

# Plot each gene
for _, row in region_sorted.iterrows():
    y = row["track"] * 0.5
    strand = row["strand"]
    color = "steelblue" if strand == "+" else "indianred"
    name = row["gene_name"][:20]  # shorten name if needed
    start = row["start"]
    end = row["end"]

    # Draw gene body with arrow
    arrow = mpatches.FancyArrow(
        start if strand == "+" else end,
        y,
        end - start if strand == "+" else start - end,
        0,
        width=0.1,
        length_includes_head=True,
        head_width=0.2,
        head_length=500,
        color=color,
        alpha=0.9
    )
    ax.add_patch(arrow)

    # Add label
    ax.text((start + end) / 2, y + 0.2, name,
            ha="center", va="bottom", fontsize=7, rotation=0)

# Layout
ax.set_xlim(args.start, args.end)
ax.set_ylim(-0.5, max(region_sorted["track"]) * 0.6 + 1)
ax.set_yticks([])
ax.set_xlabel("Genomic Position (bp)")
ax.set_title(f"Genes in {args.seqid}:{args.start}-{args.end}", fontsize=12, weight="bold")

plt.tight_layout()
plt.savefig(args.out_png)
plt.savefig(args.out_pdf)