###########################################################################
#
# ███████╗███████╗████████╗██╗   ██╗██████╗ 
# ██╔════╝██╔════╝╚══██╔══╝██║   ██║██╔══██╗
# ███████╗█████╗     ██║   ██║   ██║██████╔╝
# ╚════██║██╔══╝     ██║   ██║   ██║██╔═══╝ 
# ███████║███████╗   ██║   ╚██████╔╝██║     
# ╚══════╝╚══════╝   ╚═╝    ╚═════╝ ╚═╝     
#                                         
# This file contains rules for importing the data
# from the shared storage and clean them up for further analysis.

###########################################################################

################################################
## Rule: import_data
## Description: This rule copies raw VCF files from the shared storage to the local working directory.
################################################

rule import_data:
    """
    This rule copies raw data from the shared storage to the local working directory.
    """
    output:
        directory("data/raw")
    log:
        err = "logs/setup/import_data.err",
        out = "logs/setup/import_data.out"
    conda:
        "../envs/setup.yaml"
    resources:
        mem_mb = 2000,
        cpus_per_task = 1,
        threads = 1,
        runtime = "10m"
    shell:
        """
        mkdir -p {output}
        rsync -av --progress /nas/FAC/FBM/DEE/mrobinso/default/D2c/mbrasovi/Banyuls_Roscoff/DNAseqVCF/ {output}/ > {log.out} 2> {log.err}
        """

################################################
## Rule: subset_vcf
## Description: This rule creates a subset of a VCF file for faster testing by keeping the first 10,000 variant lines and all headers.
################################################

rule subset_vcf:
    """
    This rule creates a subset of a VCF file for faster testing by keeping the first 10,000 variant lines and all headers.
    """
    input:
        vcf = "data/raw/ShortVariants_HardCallableFiltered.{chromosomes}.vcf.gz"
    output:
        vcf = "data/subset/ShortVariants_HardCallableFiltered.{chromosomes}.vcf.gz"
    log:
        err = "logs/setup/subset_vcf_{chromosomes}.err",
        out = "logs/setup/subset_vcf_{chromosomes}.out"
    conda:
        "../envs/setup.yaml"
    resources:
        mem_mb = 4000,
        cpus_per_task = 1,
        threads = 1,
        runtime = "10m"
    shell:
        """
        (
            bcftools view -h {input.vcf} || true
            bcftools view -H {input.vcf} | head -n 10000 || true
        ) | bcftools view -Oz -o {output.vcf}

        bcftools index {output.vcf} > {log.out} 2> {log.err}
        """

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
        err = "logs/setup/copy_metadata_and_reference.err",
        out = "logs/setup/copy_metadata_and_reference.out"
    conda:
        "../envs/setup.yaml"
    resources:
        mem_mb = 2000,
        cpus_per_task = 1,
        threads = 1,
        runtime = "10m"
    shell:
        """
        cp {input.metadata} {output[0]} && \
        cp {input.reference} {output[1]} > {log.out} 2> {log.err}
        """

################################################
## Rule: copy_slrfinder_scripts
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
        err = "logs/setup/copy_slrfinder_scripts.err",
        out = "logs/setup/copy_slrfinder_scripts.out"
    resources:
        mem_mb = 1000,
        cpus_per_task = 1,
        threads = 1,
        runtime = "2m"
    shell:
        """
        cp {input.functions} {output.functions_out} && \
        cp {input.scripts} {output.scripts_out} > {log.out} 2> {log.err}
        """