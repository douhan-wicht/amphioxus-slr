#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(scales)
  library(dplyr)
})

# Command-line options
option_list <- list(
  make_option("--input", type = "character", help = "Path to .rds file"),
  make_option("--out_png", type = "character", help = "Output PNG path"),
  make_option("--out_pdf", type = "character", help = "Output PDF path"),
  make_option("--out_svg", type = "character", help = "Output SVG path")
)
opt <- parse_args(OptionParser(option_list = option_list))

# Load RDS
data_all <- readRDS(opt$input)

# Validate structure
if (!all(c("cluster", "data") %in% colnames(data_all))) {
  stop("Expected 'cluster' and 'data' columns in the input dataframe.")
}

# Filter for cluster 4471
target_row <- data_all %>% filter(cluster == 4471)
if (nrow(target_row) == 0) {
  stop("Cluster 4471 not found in 'cluster' column.")
}

# Extract nested dataframe
cluster_data <- target_row$data[[1]]

# Validate columns
required_cols <- c("PC1", "Het", "sex")
missing <- setdiff(required_cols, colnames(cluster_data))
if (length(missing) > 0) {
  stop(paste("Missing columns in nested dataframe:", paste(missing, collapse = ", ")))
}

# Rescale
cluster_data$PC1_scaled <- rescale(cluster_data$PC1)
cluster_data$Het_scaled <- rescale(cluster_data$Het)

# Create a readable sex label column for the legend
cluster_data$sex_label <- case_when(
  cluster_data$sex == "Male" ~ "Male",
  cluster_data$sex == "Female" ~ "Female",
  TRUE ~ "Unknown"
)

# Assign Seaborn colorblind palette to sex groups
sex_colors <- c(
  "Male" = "#DE8F05",
  "Female" = "#0173B2",
  "Unknown" = "#029E73"
)

# Plot
p <- ggplot(cluster_data, aes(x = PC1_scaled, y = Het_scaled, color = sex_label)) +
  geom_point(alpha = 1, size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  scale_color_manual(values = sex_colors, name = "Sex") +
  labs(
    x = "PC1",
    y = "Heterozygosity"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA)
  )

# Save with transparent background
ggsave(opt$out_png, plot = p, width = 7, height = 5, dpi = 300, bg = "transparent")
ggsave(opt$out_pdf, plot = p, width = 7, height = 5, bg = "transparent")
ggsave(opt$out_svg, plot = p, width = 7, height = 5, bg = "transparent")
