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
        mem_mb = 8000, cpus_per_task = 1, threads = 1, runtime = "15m"
    shell:
        """
        awk '$3 == "gene"' {input.gff} |
        cut -f9 |
        sed -n 's/.*ID=\([^;]*\).*/\1/p' |
        sort -u > {output.genes}
        """

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
        end = 6177987 + 50000
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

rule extract_gene_info:
    output:
        raw = "data/go/gene_info.gz",
        filtered = "data/go/gene_info_7740.tsv",
        map = "data/go/locus_tag_to_geneid.tsv"
    log:
        out = "logs/go/extract_gene_info.out",
        err = "logs/go/extract_gene_info.err"
    conda:
        "../envs/go.yaml"
    shell:
        """
        # Step 1: Download full gene_info file from NCBI
        wget -N -O {output.raw} ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz \
            > {log.out} 2> {log.err}

        # Step 2: Filter only rows for amphioxus species (TaxID: 7740) + include header
        zcat {output.raw} | awk -F'\t' 'NR==1 || $1 == "7740"' > {output.filtered}

        # Step 3: Extract locus_tag and GeneID
        awk -F '\t' 'BEGIN{{OFS="\t"}} NR>1 && $7 != "-" {{print $7, $2}}' {output.filtered} > {output.map}
        """


rule build_locus_to_geneid_map:
    input:
        gff = "data/annotation/genomic.gff",
        gene_info = "data/go/gene_info_7740.tsv"
    output:
        map = "data/annotation/locus_to_geneid_map.tsv"
    conda:
        "../envs/go.yaml"
    log:
        out = "logs/go/build_locus_to_geneid_map.out",
        err = "logs/go/build_locus_to_geneid_map.err"
    resources:
        mem_mb = 2000, cpus_per_task = 1, threads = 1, runtime = "5m"
    shell:
        """
        python workflow/scripts/go/build_locus_to_geneid_map.py \
            --gff {input.gff} \
            --gene_info {input.gene_info} \
            --output {output.map} \
            > {log.out} 2> {log.err}
        """

rule convert_target_genes:
    input:
        genes = "results/go/genes_in_region.txt",
        map = "data/annotation/locus_to_geneid_map.tsv"
    output:
        ids = "results/go/target_geneids.txt"
    log:
        out = "logs/go/convert_target_genes.out",
        err = "logs/go/convert_target_genes.err"
    conda:
        "../envs/go.yaml"
    resources:
        mem_mb = 8000, cpus_per_task = 1, threads = 1, runtime = "10m"
    shell:
        """
        python workflow/scripts/go/convert_gene_list_to_ids.py \
            --gene_list {input.genes} \
            --mapping {input.map} \
            --output {output.ids} \
            > {log.out} 2> {log.err}
        """

rule convert_background_genes:
    input:
        genes = "data/annotation/all_genes.txt",
        map = "data/annotation/locus_to_geneid_map.tsv"
    output:
        ids = "data/annotation/all_geneids.txt"
    log:
        out = "logs/go/convert_background_genes.out",
        err = "logs/go/convert_background_genes.err"
    conda:
        "../envs/go.yaml"
    resources:
        mem_mb = 8000, cpus_per_task = 1, threads = 1, runtime = "10m"
    shell:
        """
        python workflow/scripts/go/convert_gene_list_to_ids.py \
            --gene_list {input.genes} \
            --mapping {input.map} \
            --output {output.ids} \
            > {log.out} 2> {log.err}
        """

rule extract_gene2go:
    output:
        raw = "data/go/gene2go.gz",
        filtered = "data/go/gene2go_7740.tsv",
        obo = "data/go/go-basic.obo"
    log:
        out = "logs/go/extract_gene2go.out",
        err = "logs/go/extract_gene2go.err"
    conda:
        "../envs/go.yaml"
    resources:
        mem_mb = 4000, cpus_per_task = 1, threads = 1, runtime = "10m"
    shell:
        """
        mkdir -p data/go logs/go

        wget -N -O {output.raw} ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz \
            >> {log.out} 2>> {log.err}

        zgrep -E "^GeneID|^7740" {output.raw} > {output.filtered}

        wget -N -O {output.obo} http://purl.obolibrary.org/obo/go/go-basic.obo \
            >> {log.out} 2>> {log.err}
        """

rule go_enrichment:
    """
    Run GO enrichment using goatools Fisher’s exact test.
    """
    input:
        gene2go = "data/go/gene2go_7740.tsv",
        obo = "data/go/go-basic.obo",
        target_genes = "results/go/target_geneids.txt",
        background_genes = "data/annotation/all_geneids.txt"
    output:
        enriched = "results/go/go_enrichment.tsv"
    log:
        out = "logs/go/enrichment.out",
        err = "logs/go/enrichment.err"
    conda:
        "../envs/go.yaml"
    resources:
        mem_mb = 8000, cpus_per_task = 1, threads = 1, runtime = "10m"
    shell:
        """
        python workflow/scripts/go/go_enrichment.py \
            --gene2go {input.gene2go} \
            --obo {input.obo} \
            --target {input.target_genes} \
            --background {input.background_genes} \
            --out {output.enriched} \
            > {log.out} 2> {log.err}
        """
