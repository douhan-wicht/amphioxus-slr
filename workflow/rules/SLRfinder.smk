###########################################################################
#
# ███████╗██╗     ██████╗ ███████╗██╗███╗   ██╗██████╗ ███████╗██████╗ 
# ██╔════╝██║     ██╔══██╗██╔════╝██║████╗  ██║██╔══██╗██╔════╝██╔══██╗
# ███████╗██║     ██████╔╝█████╗  ██║██╔██╗ ██║██║  ██║█████╗  ██████╔╝
# ╚════██║██║     ██╔══██╗██╔══╝  ██║██║╚██╗██║██║  ██║██╔══╝  ██╔══██╗
# ███████║███████╗██║  ██║██║     ██║██║ ╚████║██████╔╝███████╗██║  ██║
# ╚══════╝╚══════╝╚═╝  ╚═╝╚═╝     ╚═╝╚═╝  ╚═══╝╚═════╝ ╚══════╝╚═╝  ╚═╝    
#                                         
# This file contains rules for running SLRfinder on the formated data.

###########################################################################

################################################
## Rule: vcf_filtering_ld_estimation
## Description: This rule filters the VCF files and estimates LD using vcftools.
## Look for PASS flag -> the ones that have passed the Marina check. (gatk calling)
################################################

rule vcf_filtering_ld_estimation:
    input:
        # The VCF files for each chromosome/contig and the reference list
        vcf="data/subset/ShortVariants_HardCallableFiltered.{chromosomes}.vcf.gz",
        reference="tmp/amphioxus/reference.list"
    output:
        # Output filtered VCF and LD edge list
        filtered_vcf="tmp/amphioxus/a15m75/amphioxus_{chromosomes}_a15m75.recode.vcf",
        ld_file="tmp/amphioxus/GenoLD.snp100/amphioxus_{chromosomes}_a15m75.geno.ld"
    log:
        err = "logs/SLRfinder/vcf_filtering_ld_estimation_{chromosomes}.err",
        out = "logs/SLRfinder/vcf_filtering_ld_estimation_{chromosomes}.out"
    conda:
        '../envs/SLRfinder.yaml'
    params:
        min_ac = 1,
        min_gq = 20,
        min_q = 30,
        maf = 0.15,
        max_missing = 0.75,
        ld_window = 100
    resources:
        mem_mb = 1000,
        cpus_per_task = 1,
        threads = 1,
        runtime = "15m"
    shell:
        """
        # Create necessary directories inside the 'tmp/amphioxus' folder if they don't exist
        mkdir -p tmp/amphioxus/a15m75
        mkdir -p tmp/amphioxus/GenoLD.snp100

        # Step 1: SNP filtering using bcftools and vcftools
        bcftools view -m2 -M2 -v snps --min-ac={params.min_ac} {input.vcf} \
        | vcftools --vcf - --minGQ {params.min_gq} --minQ {params.min_q} --maf {params.maf} --max-missing {params.max_missing} \
        --recode --recode-INFO-all --out tmp/amphioxus/a15m75/amphioxus_{wildcards.chromosomes}_a15m75 \
        > {log.out} 2> {log.err}

        # Step 2: LD estimation using vcftools
        vcftools --vcf tmp/amphioxus/a15m75/amphioxus_{wildcards.chromosomes}_a15m75.recode.vcf --geno-r2 --ld-window {params.ld_window} \
        --out tmp/amphioxus/GenoLD.snp100/amphioxus_{wildcards.chromosomes}_a15m75 \
        > {log.out} 2> {log.err}
        """

################################################
## Rule: SLRfinder_main
## Description: This rule runs the SLRfinder analysis on the filtered VCF files.
################################################

rule SLRfinder_main:
    input:
        "tmp/amphioxus/SLRfinder_scripts.R",
        "tmp/amphioxus/SLRfinder_functions.r"
    output:
        directory("tmp/amphioxus/LD8.5cl20")
    log:
        err = "logs/SLRfinder/SLRfinder_main.err",
        out = "logs/SLRfinder/SLRfinder_main.out"
    conda:
        '../envs/SLRfinder.yaml'
    resources:
        mem_mb = 16000,
        cpus_per_task = 1,
        threads = 1,
        runtime = "1h"
    shell:
        """
        cd tmp/amphioxus
        Rscript SLRfinder_scripts.R > ./../../{log.out} 2> ./../../{log.err}
        """
