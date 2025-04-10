library(data.table)
library(igraph)
library(ggpubr)
library(SNPRelate)

mydata = "amphioxus"
min_LD = 0.85
min.cl.size = 20
myranks = c("Dext_var_rank", "R2_rank", "nSNPs_rank", "chi2_rank")
sex_info = TRUE
sex_filter = 0.1
my_sex_ratio = c(0.5, 0.5)
ncores = 1

print("Reading data information...")

sif = read.csv(paste0(mydata, ".csv"))
LG = read.table("reference.list", header = FALSE)
names(LG) = c("chr", "lg")

source("SLRfinder_functions.r")
print("Sourced SLR functions.")

# Step 1: Get LD clusters
print("Step 1: Getting LD clusters")

dir.create(paste0("LD", min_LD*10, "cl", min.cl.size), showWarnings = FALSE)
setwd(paste0("LD", min_LD*10, "cl", min.cl.size))
dir.create("whitelist", showWarnings = FALSE)

print("Loading LD data...")

data_cls <- NULL

for (i in 1:nrow(LG)) {
  chr = LG[i, "chr"]
  lg = LG[i, "lg"]

  print(paste0("Processing chromosome ", chr, " (", lg, ")..."))

  ld_file = paste0("../GenoLD.snp100/", mydata, "_", lg, "_a15m75.geno.ld")
  if (!file.exists(ld_file)) {
    cat("⚠️  Skipping", chr, "- missing LD file\n")
    next
  }

  data = read.table(ld_file, header = TRUE)
  names(data) = c("CHR", "from", "to", "N_INDV", "r2")

  out = get_single_LD_cluster(data, min_LD = min_LD, min.cl.size = min.cl.size)

  if (!is.null(out) && nrow(out) > 0) {
    position = as.data.frame(unlist(out$SNPs))
    position = cbind(rep(chr, sum(out$nSNPs)), position)
    whitelist_path = paste0("whitelist/position.", lg, ".list")
    write.table(position, whitelist_path, sep = "\t", quote = FALSE, row.names = FALSE)
    data_cls <- rbind(data_cls, out)
    cat("✔️  Wrote whitelist for", chr, "\n")
  } else {
    cat("⚠️  No LD clusters found for", chr, "\n")
  }
}

if (is.null(data_cls)) {
  stop("❌ No LD clusters found. Exiting.")
}

data_cls$SNPs = apply(data_cls, 1, function(cl) {
  paste0(cl$chr, "_", cl$SNPs)
})

print(paste0("Total number of LD clusters: ", nrow(data_cls)))
saveRDS(data_cls, file = "data_cls.rds")

# Step 2: Generate 012 matrices
dir.create("file012", showWarnings = FALSE)

print("Generating 012 matrices...")

for (i in 1:nrow(LG)) {
  lg = LG[i, "lg"]
  whitelist_file = paste0("whitelist/position.", lg, ".list")
  vcf_file = paste0("../a15m75/", mydata, "_", lg, "_a15m75.recode.vcf")
  out_file = paste0("file012/", mydata, "_", lg, "_a15m75_LD", min_LD, "cl", min.cl.size)

  print(paste0("Processing chromosome ", lg, "..."))

  if (file.exists(whitelist_file)) {
    cmd = paste(
      "vcftools",
      "--vcf", vcf_file,
      "--positions", whitelist_file,
      "--012",
      "--out", out_file
    )
    system(cmd)
    cat("✔️  Generated 012 files for", lg, "\n")
  } else {
    cat("⚠️  Skipping", lg, "- whitelist file not found\n")
  }
}

# Step 3: Load genotypes and map
print("Step 2: Processing LD clusters")

data_cls = readRDS("data_cls.rds")

files <- list.files("file012", full.names = TRUE)
indv_files <- files[grep(".indv", files)]
pos_files <- files[grep(".pos", files)]
GT_files <- files[!grepl(".log", files) & !grepl(".indv", files) & !grepl(".pos", files)]

map <- rbindlist(lapply(pos_files, function(pos_file) {
  pos <- fread(pos_file, sep = "\t")
  if (ncol(pos) > 1) {
    colnames(pos) <- c("Chr", "Pos")
    pos$SNP <- paste0(pos$Chr, "_", pos$Pos)
  }
  if (ncol(pos) == 0) pos = NULL
  return(pos)
}))

GT <- do.call(cbind, lapply(GT_files, function(gt_file) {
  gt.matrix = as.matrix(fread(gt_file)[, -1])
  if (nrow(gt.matrix) == 0) return(NULL)
  return(gt.matrix)
}))
GT[GT == -1] <- NA

indv <- fread(indv_files[1], header = FALSE)
pop_info <- sif[order(factor(sif$SampleID, levels = indv$V1)), ]
if (!all(indv$V1 == pop_info$SampleID)) stop("❌ Individual order mismatch.")

ind <- pop_info$SampleID
pop <- pop_info$Population

save(data_cls, GT, map, ind, pop, file = "GT.RData")

data_all = get_data_output(data_cls, GT, map, pop, sex_info, heterog_homog = my_sex_ratio, cores = ncores)
saveRDS(data_all, "data_all.rds")

# Step 4: Identify SLR candidates
print("Step 3: Identify SLR candidates")

if (sex_info) {
  print(paste0("Filtering clusters by sex (≤ ", sex_filter*100, "% misgrouped)"))
  data_sex = data_all[data_all$Sex_g <= sex_filter, ]

  if (nrow(data_sex) > 0) {
    myindex = length(grep("_", data_sex[1, "chr"])) + 2
    data_sex[, region := apply(data_sex, 1, function(x) {
      paste0(x$chr, ":", paste(range(as.numeric(do.call(rbind, strsplit(x$SNPs, "_", fixed = TRUE))[, myindex])), collapse = "-"))
    })]


    pdf(paste0(mydata, "_sexg.pdf"))
    for (i in 1:nrow(data_sex)) {
      d = as.data.frame(data_sex$data[i])
      region = data_sex$region[i]
      title = paste0(region, "\nnSNPs=", data_sex$nSNPs[i])

      print(ggplot(d, aes(x = PC_scaled, y = Het)) +
              geom_point(aes(color = sex), alpha = 0.6, size = 2.5) +
              geom_smooth(method = "lm", se = FALSE, col = "black") +
              theme_bw() +
              labs(title = title))
    }
    dev.off()

    non_list_cols <- names(which(sapply(data_sex, function(col) !is.list(col))))
    data_sex_export <- data_sex[, ..non_list_cols]
    write.csv(data_sex_export, "sex_filter.csv", row.names = FALSE)
  } else {
    print(paste0("No cluster passed sex_g ≤ ", sex_filter))
  }
}

print("Identifying candidates by rank...")

cand_regions <- get_candidate_regions(data_all, ranks = myranks, nPerm = 10000, cores = ncores)
saveRDS(cand_regions, "cand_regions.rds")

# Final visualization
list2env(cand_regions, globalenv())
alpha = 0.05
lambda <- lm(obs ~ exp + 0, qq_data)$coefficients
qq_data$col <- ifelse(data_out$p_gc_adj < alpha, "#ff9727", "grey40")

pdf(paste0(mydata, "_LD", min_LD, "cl", min.cl.size, ".pdf"), width = 6, height = 4)
print(ggplot(qq_data, aes(x = exp, y = obs)) +
        geom_point(col = qq_data$col) +
        theme_bw() +
        labs(title = paste0("min_LD=", min_LD, ", min.cl.size=", min.cl.size, "\nλ=", round(lambda, 2)),
             x = "Expected -log10(P)", y = "Observed -log10(P)") +
        geom_abline(slope = 1, intercept = 0, linewidth = 0.5) +
        geom_smooth(method = "lm", col = "#ff9727", linewidth = 0.5))

# Plot candidate regions

print("Plotting candidate regions...")

for (r in unique(PCA_het_data$region)) {
  pca = PCA_het_data[PCA_het_data$region == r, ]
  label = strsplit(unique(pca$label), " ")[[1]]
  chr = strsplit(label[1], ":")[[1]][1]
  lg = LG[LG$chr == chr, "lg"]
  title = paste0(sub(chr, lg, label[1]), "\n", label[2], " ", label[3])

  if (sex_info) {
    print(ggplot(pca, aes(PC_scaled, Het)) +
            geom_point(aes(color = sex), alpha = 0.6, size = 2.5) +
            geom_smooth(method = "lm", se = FALSE, col = "black") +
            theme_bw() + labs(title = title))
  } else {
    print(ggplot(pca, aes(PC_scaled, Het)) +
            geom_point(aes(color = Pop), alpha = 0.6, size = 2.5) +
            geom_smooth(method = "lm", se = FALSE, col = "black") +
            theme_bw() + labs(title = title))
  }
}
dev.off()

setwd("../../")
print("🎉 SLRfinder pipeline completed.")
