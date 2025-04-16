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

################################################
## Rule: extract_background_gene_list
## Description: Extract background gene list from GFF annotation file.
################################################

rule extract_background_gene_list:
    """
    Extract background gene list from GFF annotation file.
    """
    input:
        gff = "data/annotation/genomic.gff"
    output:
        genes = "data/annotation/all_genes.txt"
    conda:
        "../envs/go.yaml"
    resources:
        mem_mb = 8000,
        cpus_per_task = 1,
        threads = 1,
        runtime = "15m"
    shell:
        """
        awk '$3 == "gene"' {input.gff} |
        cut -f9 |
        sed -n 's/.*ID=\\([^;]*\\).*/\\1/p' |
        sort -u > {output.genes}
        """

################################################
## Rule: extract_gene_list
## Description: Extract gene names from GFF in a region of interest.
################################################

rule extract_gene_list:
    """
    Extract gene names from GFF in a region of interest.
    """
    input:
        gff = "data/annotation/genomic.gff"
    output:
        gene_list = "results/go/genes_in_region.txt"
    log:
        out = "logs/go/extract_genes.out",
        err = "logs/go/extract_genes.err"
    conda:
        "../envs/go.yaml"
    params:
        seqid = "OV696689.1",
        start = 6132346 - 50000,
        end = 6174195 + 50000
    shell:
        """
        python workflow/scripts/go/extract_gene_list.py \
            --gff {input.gff} \
            --seqid {params.seqid} \
            --start {params.start} \
            --end {params.end} \
            --output {output.gene_list} \
            > {log.out} 2> {log.err}
        """


