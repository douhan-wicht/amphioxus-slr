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
