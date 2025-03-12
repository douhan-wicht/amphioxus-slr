name: SLRfinder
channels:
    - conda-forge
    - bioconda
    - defaults
dependencies:
    - r-base=4.2.3
    - r-igraph
    - r-data.table
    - r-ggplot2
    - r-ggpubr
    - r-cowplot
    - r-parallel
    - r-BiocManager  # This ensures that BiocManager is available for installing Bioconductor packages
    - r-SNPRelate  # This will work if the Bioconductor package is available in conda
    - r-snpStats  # If you're working with SNP-related data, this could be another useful package
post-activate:
    - Rscript -e "if (!require('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install('SNPRelate')"
