import argparse
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
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
female_color = palette[0]
male_color = palette[3]
diff_color = palette[2]
region_color = palette[6]
arrow_color = palette[0]

# -----------------------------
# LEGEND ELEMENTS
# -----------------------------
legend_elements = [
    mpatches.Patch(color=female_color, label="Females"),
    mpatches.Patch(color=male_color, label="Males"),
    mpatches.Patch(color=diff_color, label="Females - Males"),
    mlines.Line2D([0], [0], color=arrow_color, lw=2,
                  marker='>', markersize=10, label="Gene Arrow", linestyle='-'),
    mpatches.Patch(color=region_color, label="Region of Interest")
]

# -----------------------------
# LEGEND-ONLY FIGURE (Vertical)
# -----------------------------
fig, ax = plt.subplots(figsize=(3.5, 3))
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
