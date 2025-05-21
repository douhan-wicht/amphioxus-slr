###########################################################################
#
#  ██████╗ ███████╗███╗   ██╗███████╗     ██████╗ ███╗   ██╗████████╗ ██████╗ ██╗      ██████╗  ██████╗██╗   ██╗
# ██╔════╝ ██╔════╝████╗  ██║██╔════╝    ██╔═══██╗████╗  ██║╚══██╔══╝██╔═══██╗██║     ██╔═══██╗██╔════╝╚██╗ ██╔╝
# ██║  ███╗█████╗  ██╔██╗ ██║█████╗      ██║   ██║██╔██╗ ██║   ██║   ██║   ██║██║     ██║   ██║██║  ███╗╚████╔╝ 
# ██║   ██║██╔══╝  ██║╚██╗██║██╔══╝      ██║   ██║██║╚██╗██║   ██║   ██║   ██║██║     ██║   ██║██║   ██║ ╚██╔╝  
# ╚██████╔╝███████╗██║ ╚████║███████╗    ╚██████╔╝██║ ╚████║   ██║   ╚██████╔╝███████╗╚██████╔╝╚██████╔╝  ██║   
#  ╚═════╝ ╚══════╝╚═╝  ╚═══╝╚══════╝     ╚═════╝ ╚═╝  ╚═══╝   ╚═╝    ╚═════╝ ╚══════╝ ╚═════╝  ╚═════╝   ╚═╝      
#                                         
# This file contains rules for doing a Gene Ontology (GO) enrichment analysis

###########################################################################

rule check_haplotype_pattern:
    input:
        table="tmp/amphioxus/amphioxus_chr4.tab"
    output:
        table="results/misc/haplotype_check_{start}_{end}.tsv",
        summary="results/misc/haplotype_check_{start}_{end}_summary.txt"
    log:
        err="logs/misc/haplotype_check_{start}_{end}.err"
    conda:
        "../envs/misc.yaml"
    params:
        start=lambda wildcards: config["region"]["start"],
        end=lambda wildcards: config["region"]["end"]
    shell:
        """
        python workflow/scripts/misc/check_haplotype_pattern.py \
            --input {input.table} \
            --start {params.start} \
            --end {params.end} \
            --output-prefix results/misc/haplotype_check_{params.start}_{params.end} \
            2> {log.err}
        """

rule check_haplotype_pattern_combined:
    input:
        table="tmp/amphioxus/amphioxus_chr4.tab",
        pos_list="results/misc/ld_cluster_snps.txt"
    output:
        table="results/misc/haplotype_check_combined.tsv",
        summary="results/misc/haplotype_check_combined_summary.txt"
    log:
        err="logs/misc/haplotype_check_combined.err"
    conda:
        "../envs/misc.yaml"
    params:
        start=lambda wildcards: config["region"]["start"],
        end=lambda wildcards: config["region"]["end"]
    shell:
        """
        python workflow/scripts/misc/check_haplotype_pattern.py \
            --input {input.table} \
            --start {params.start} \
            --end {params.end} \
            --pos-list {input.pos_list} \
            --output-prefix results/misc/haplotype_check_combined \
            2> {log.err}
        """

rule normalize_ld_clusters:
    input:
        clusters="tmp/amphioxus/LD8.5cl20/candidates.csv"
    output:
        norm_table="results/misc/normalized_ld_clusters.tsv"
    log:
        err="logs/misc/normalize_ld_clusters.err"
    resources:
        mem_mb = 8000,
        cpus_per_task = 1,
        threads = 1,
        runtime = "30m"
    conda:
        "../envs/misc.yaml"
    shell:
        """
        python workflow/scripts/misc/normalize_clusters_by_length.py \
            --input {input.clusters} \
            --output {output.norm_table} \
            2> {log.err}
        """


# vcf="data/raw/ShortVariants_HardCallableFiltered.{chromosome}.vcf.gz"
rule count_raw_snps:
    input:
        vcf="tmp/amphioxus/a15m75/amphioxus_{chromosome}_a15m75.recode.vcf"
    output:
        snp_count="results/misc/snp_counts/ShortVariants_HardCallableFiltered.{chromosome}.txt"
    log:
        err="logs/misc/snp_counts_{chromosome}.err"

    resources:
        mem_mb = 10000,
        cpus_per_task = 1,
        threads = 1,
        runtime = "10m"
    conda:
        "../envs/misc.yaml"
    shell:
        """
        mkdir -p results/misc/snp_counts

        vcftools \
            --gzvcf {input.vcf} \
            --remove-filtered-all \
            --remove-indels \
            --recode \
            --stdout \
        | grep -v '^#' \
        | wc -l \
        > {output.snp_count} 2> {log.err}
        """

rule merge_snp_counts:
    input:
        expand("results/misc/snp_counts/ShortVariants_HardCallableFiltered.{chromosome}.txt", chromosome=config["CHROMOSOMES"])
    output:
        summary="results/misc/snp_counts/summary_snp_counts.tsv"
    log:
        err="logs/misc/merge_snp_counts.err"
    conda:
        "../envs/misc.yaml"
    shell:
        """
        python workflow/scripts/misc/merge_snp_counts.py \
            --input {input} \
            --output {output.summary} \
            2> {log.err}
        """