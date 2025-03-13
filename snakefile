###########################################################################

#  ███████╗███╗   ██╗ █████╗ ██╗  ██╗███████╗███████╗██╗██╗     ███████╗
#  ██╔════╝████╗  ██║██╔══██╗██║ ██╔╝██╔════╝██╔════╝██║██║     ██╔════╝
#  ███████╗██╔██╗ ██║███████║█████╔╝ █████╗  █████╗  ██║██║     █████╗  
#  ╚════██║██║╚██╗██║██╔══██║██╔═██╗ ██╔══╝  ██╔══╝  ██║██║     ██╔══╝  
#  ███████║██║ ╚████║██║  ██║██║  ██╗███████╗██║     ██║███████╗███████╗
#  ╚══════╝╚═╝  ╚═══╝╚═╝  ╚═╝╚═╝  ╚═╝╚══════╝╚═╝     ╚═╝╚══════╝╚══════╝
#
# This file is used to produce all the output files necessary to reproduce the paper's results.
# Workflow:
#   1. Import the project wide configuration for snakemake.
#   2. Include local rules to import the data set (will have to be changed to enhance reproducibility).
#   3. Include rules to process and analyse the data.
#   4. Include the master rule to produce the final output files.

###########################################################################

# Project wide configuration
configfile: "config/amphioxus-slr.yaml"

# Include local rules
localrules: import_data

# Include the other rules
include: "rules/setup.smk",
include: "rules/SLRfinder.smk"

# Include a master rule to produce all the final output files
rule all:
    input:
        "scripts/SLRfinder/amphioxus/amphioxus.csv", # Metadata
        "scripts/SLRfinder/amphioxus/reference.list", # Reference list
        expand("scripts/SLRfinder/amphioxus/a15m75/amphioxus_{chromosomes}_a15m75.recode.vcf", chromosomes=config["CHROMOSOMES"]),
        expand("scripts/SLRfinder/amphioxus/GenoLD.snp100/amphioxus_{chromosomes}_a15m75.geno.ld", chromosomes=config["CHROMOSOMES"])