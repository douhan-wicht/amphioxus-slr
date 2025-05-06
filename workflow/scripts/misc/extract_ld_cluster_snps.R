args <- commandArgs(trailingOnly = TRUE)
rds_path <- args[1]
output_path <- args[2]

# Load RDS file
data_all <- readRDS(rds_path)

# Check structure
if (!"cluster" %in% names(data_all) || !"SNPs" %in% names(data_all)) {
  stop("Expected columns 'cluster' and 'SNPs' not found.")
}

# Get SNPs from cluster 4471
snps_raw <- data_all$SNPs[data_all$cluster == 4471][[1]]

# Extract numeric POS from e.g., "chr4_6142346"
pos <- as.integer(sub(".*_", "", snps_raw))

# Save to output file
write.table(sort(unique(pos)), file=output_path, quote=FALSE, row.names=FALSE, col.names=FALSE)
