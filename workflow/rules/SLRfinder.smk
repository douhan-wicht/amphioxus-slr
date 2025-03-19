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
## Rule: vcf_filtering_ld_estimation
## Description: This rule filters the VCF files and estimates LD using vcftools.
################################################

# rule vcf_filtering_ld_estimation:
#     input:
#         # The VCF files for each chromosome/contig and the reference list
#         vcf="data/DNAseqVCF/ShortVariants_HardCallableFiltered.{chromosomes}.vcf.gz",
#         reference="scripts/SLRfinder/amphioxus/reference.list"
#     output:
#         # Output filtered VCF and LD edge list
#         filtered_vcf="scripts/SLRfinder/amphioxus/a15m75/amphioxus_{chromosomes}_a15m75.recode.vcf",
#         ld_file="scripts/SLRfinder/amphioxus/GenoLD.snp100/amphioxus_{chromosomes}_a15m75.geno.ld"
#     log:
#         err = "logs/SLRfinder/vcf_filtering_ld_estimation_{chromosomes}.err",
#         out = "logs/SLRfinder/vcf_filtering_ld_estimation_{chromosomes}.out"
#     conda:
#         '../envs/SLRfinder.yaml'
#     params:
#         time = "00:30:00",
#         name = "copy_metadata_and_reference",
#         threads = 20,
#         mem = "8G",
#         # Additional parameters
#         min_ac = 1,
#         min_gq = 20,
#         min_q = 30,
#         maf = 0.15,
#         max_missing = 0.75,
#         ld_window = 100
#     shell:
#         """
#         # Create necessary directories inside the 'scripts/SLRfinder/amphioxus' folder if they don't exist
#         mkdir -p scripts/SLRfinder/amphioxus/a15m75
#         mkdir -p scripts/SLRfinder/amphioxus/GenoLD.snp100

#         # Get chromosome name and LG from the reference list using the SLURM_ARRAY_TASK_ID
#         chr=$(sed -n ${{SLURM_ARRAY_TASK_ID}}p {input.reference} | awk '{{print $1}}')
#         lg=$(sed -n ${{SLURM_ARRAY_TASK_ID}}p {input.reference} | awk '{{print $2}}')
#         ind=amphioxus

#         # Step 1: SNP filtering using bcftools and vcftools
#         bcftools view -m2 -M2 -v snps --min-ac={params.min_ac} {input.vcf} \
#         | vcftools --vcf - --minGQ {params.min_gq} --minQ {params.min_q} --maf {params.maf} --max-missing {params.max_missing} \
#         --recode --recode-INFO-all --out {output.filtered_vcf} \
#         > {log.out} 2>> {log.err}

#         # Step 2: LD estimation using vcftools
#         vcftools --vcf {output.filtered_vcf} --geno-r2 --ld-window {params.ld_window} \
#         --out {output.ld_file} \
#         > {log.out} 2>> {log.err}
#         """
