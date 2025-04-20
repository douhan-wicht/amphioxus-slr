import argparse
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

# -----------------------------
# ARGPARSE
# -----------------------------
parser = argparse.ArgumentParser(description="Generate vertical legend for heterozygosity + gene track plots.")
parser.add_argument("--out_png", required=True, help="Output PNG path for legend")
parser.add_argument("--out_pdf", required=True, help="Output PDF path for legend")
parser.add_argument("--out_svg", required=True, help="Output SVG path for legend")
args = parser.parse_args()

# -----------------------------
# COLORS
# -----------------------------
sns.set_palette("colorblind")
palette = sns.color_palette("colorblind")
female_color = palette[0]   # blue
male_color = palette[3]     # orange/green
diff_color = palette[2]     # purple
region_color = palette[6]   # pink/red
arrow_color = palette[0]    # blue

# -----------------------------
# LEGEND ELEMENTS
# -----------------------------
legend_elements = [
    mpatches.Patch(color=female_color, label="Females"),
    mpatches.Patch(color=male_color, label="Males"),
    mpatches.Patch(color=diff_color, label="Females - Males"),
    mpatches.Patch(color=arrow_color, label="Gene Arrow"),
    mpatches.Patch(color=region_color, label="Region of Interest")
]

# -----------------------------
# LEGEND-ONLY FIGURE (Vertical)
# -----------------------------
fig, ax = plt.subplots(figsize=(3.5, 3))  # Adjust size if needed
ax.axis("off")

legend = ax.legend(
    handles=legend_elements,
    loc="center left",
    frameon=False,
    fontsize=10
)

# Save all formats
plt.savefig(args.out_png, bbox_inches="tight")
plt.savefig(args.out_pdf, bbox_inches="tight")
plt.savefig(args.out_svg, bbox_inches="tight")
