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
## Rule: copy_metadata_and_reference
## Description: This rule copies the metadata and reference list to the SLRfinder folder.
################################################

rule copy_metadata_and_reference:
    input:
        metadata="data/metadata/metadata.csv",
        reference="data/metadata/reference.list"
    output:
        "tmp/amphioxus/amphioxus.csv",
        "tmp/amphioxus/reference.list"
    log:
        err = "logs/SLRfinder/copy_metadata_and_reference.err",
        out = "logs/SLRfinder/copy_metadata_and_reference.out"
    conda:
        "../envs/SLRfinder.yaml"
    resources:
        mem_mb = 2000,
        cpus_per_task = 1,
        threads = 1,
        runtime = "10m"
    shell:
        "cp {input.metadata} {output[0]} && cp {input.reference} {output[1]} >> {log.out} 2>> {log.err}"

################################################
## Rule: copy_slrfinder
## Description: This rule copies the SLRfinder functions to the tmp folder.
################################################

rule copy_slrfinder_scripts:
    input:
        functions = "workflow/scripts/SLRfinder/SLRfinder_functions.r",
        scripts = "workflow/scripts/SLRfinder/SLRfinder_scripts.R"
    output:
        functions_out = "tmp/amphioxus/SLRfinder_functions.r",
        scripts_out = "tmp/amphioxus/SLRfinder_scripts.R"
    log:
        err = "logs/SLRfinder/copy_slrfinder_scripts.err",
        out = "logs/SLRfinder/copy_slrfinder_scripts.out"
    conda:
        "../envs/SLRfinder.yaml"
    resources:
        mem_mb = 1000,
        cpus_per_task = 1,
        threads = 1,
        runtime = "5m"
    shell:
        """
        mkdir -p tmp/amphioxus && \
        cp {input.functions} {output.functions_out} && \
        cp {input.scripts} {output.scripts_out}
        """

################################################
## Rule: vcf_filtering_ld_estimation
## Description: This rule filters the VCF files and estimates LD using vcftools.
################################################

rule vcf_filtering_ld_estimation:
    input:
        # The VCF files for each chromosome/contig and the reference list
        vcf="data/raw/ShortVariants_HardCallableFiltered.{chromosomes}.vcf.gz",
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
