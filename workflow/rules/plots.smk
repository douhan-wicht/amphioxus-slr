###########################################################################
#
# ██████╗ ██╗      ██████╗ ████████╗███████╗
# ██╔══██╗██║     ██╔═══██╗╚══██╔══╝██╔════╝
# ██████╔╝██║     ██║   ██║   ██║   ███████╗
# ██╔═══╝ ██║     ██║   ██║   ██║   ╚════██║
# ██║     ███████╗╚██████╔╝   ██║   ███████║
# ╚═╝     ╚══════╝ ╚═════╝    ╚═╝   ╚══════╝    
#                                         
# This file contains rules for plotting the results.

###########################################################################

################################################
## Rule: karyotype
## Description: This rule plots a karyotype of the genome and colours our regions of interest
## (LD clusters) based on the percentage of missgrouped individuals.
################################################

rule karyotype:
    """
    This rule plots a karyotype of the genome and colours our regions of interest (LD clusters)
    based on the percentage of missgrouped individuals.
    """
    input:
        candidates = "tmp/amphioxus/LD8.5cl20/candidates.csv",
        sex_filter = "tmp/amphioxus/LD8.5cl20/sex_filter.csv"
    output:
        png = "results/plots/karyotype.png",
        pdf = "results/plots/karyotype.pdf",
        svg = "results/plots/karyotype.svg"
    log:
        err = "logs/plots/karyotype.err",
        out = "logs/plots/karyotype.out"
    conda:
        "../envs/plots.yaml"
    resources:
        mem_mb = 2000,
        cpus_per_task = 1,
        threads = 1,
        runtime = "2m"
    shell:
        """
        python workflow/scripts/plots/karyotype.py \
            --candidates {input.candidates} \
            --sex_filter {input.sex_filter} \
            --out_png {output.png} \
            --out_pdf {output.pdf} \
            --out_svg {output.svg} \
            > {log.out} 2> {log.err}
        """

################################################
## Rule: gene_region_plot
## Description: This rule plots genes in a specified genomic region from a GFF3 file.
################################################

rule gene_region_plot:
    """
    Plot genes in a specified genomic region from a GFF3 file.
    """
    input:
        gff = "data/annotation/genomic.gff"
    output:
        png = "results/plots/gene_region.png",
        pdf = "results/plots/gene_region.pdf",
        svg = "results/plots/gene_region.svg"
    params:
        seqid = "OV696689.1",
        start = 6142346,
        end = 6177987
    log:
        out = "logs/plots/gene_region.out",
        err = "logs/plots/gene_region.err"
    conda:
        "../envs/plots.yaml"
    resources:
        mem_mb = 2000,
        cpus_per_task = 1,
        threads = 1,
        runtime = "2m"
    shell:
        """
        python workflow/scripts/plots/gene_region_plot.py \
            --gff {input.gff} \
            --seqid {params.seqid} \
            --start {params.start} \
            --end {params.end} \
            --out_png {output.png} \
            --out_pdf {output.pdf} \
            --out_svg {output.svg} \
            > {log.out} 2> {log.err}
        """

################################################
## Rule: vcf_to_tab
## Description: Convert VCF to tabular format using GATK VariantsToTable.
################################################

rule vcf_to_tab:
    """
    Convert VCF to tabular format using GATK VariantsToTable for each chromosome.
    """
    input:
        vcf = "tmp/amphioxus/a15m75/amphioxus_{chrom}_a15m75.recode.vcf"
    output:
        tab = "tmp/amphioxus/amphioxus_{chrom}.tab"
    log:
        out = "logs/plots/vcf_to_tab_{chrom}.out",
        err = "logs/plots/vcf_to_tab_{chrom}.err"
    conda:
        "../envs/plots.yaml"
    resources:
        mem_mb = 8000,
        cpus_per_task = 1,
        threads = 1,
        runtime = "20m"
    shell:
        """
        gatk VariantsToTable \
            -V {input.vcf} \
            -O {output.tab} \
            -F CHROM -F POS -F REF -F ALT \
            -GF GT \
            > {log.out} 2> {log.err}
        """

rule concat_tabs:
    input:
        tabs = expand("tmp/amphioxus/amphioxus_{chrom}.tab", chrom=config["CHROMOSOMES"])
    output:
        combined = "tmp/amphioxus/amphioxus_all.tab"
    log:
        out = "logs/plots/concat_tabs.out",
        err = "logs/plots/concat_tabs.err"
    shell:
        """
        head -n 1 {input.tabs[0]} > {output.combined}
        tail -n +2 -q {input.tabs} >> {output.combined}
        > {log.out} 2> {log.err}
        """

################################################
## Rule: heterozygosity_plot
## Description: Plot smoothed heterozygosity by sex from a VCF tab file.
################################################

rule heterozygosity_plot:
    """
    Plot smoothed heterozygosity by sex from a VCF tab file.
    """
    input:
        tab = "tmp/amphioxus/amphioxus_chr4.tab"
    output:
        png = "results/plots/heterozygosity_plot.png",
        pdf = "results/plots/heterozygosity_plot.pdf",
        svg = "results/plots/heterozygosity_plot.svg"
    log:
        out = "logs/plots/heterozygosity_plot.out",
        err = "logs/plots/heterozygosity_plot.err"
    conda:
        "../envs/plots.yaml"
    params:
        region_start = 6142346,
        region_end = 6177987
    shell:
        """
        python workflow/scripts/plots/heterozygosity_plot.py \
            --input {input.tab} \
            --out_png {output.png} \
            --out_pdf {output.pdf} \
            --out_svg {output.svg} \
            --region_start {params.region_start} \
            --region_end {params.region_end} \
            > {log.out} 2> {log.err}
        """

################################################
## Rule: heterozygosity_raw_plot
## Description: Plot raw (non-smoothed) heterozygosity by sex using scatter points.
################################################

rule heterozygosity_raw_plot:
    """
    Plot raw (non-smoothed) heterozygosity by sex using scatter points.
    """
    input:
        tab = "tmp/amphioxus/amphioxus_chr4.tab"
    output:
        png = "results/plots/heterozygosity_raw.png",
        pdf = "results/plots/heterozygosity_raw.pdf",
        svg = "results/plots/heterozygosity_raw.svg"
    log:
        out = "logs/plots/heterozygosity_raw_plot.out",
        err = "logs/plots/heterozygosity_raw_plot.err"
    conda:
        "../envs/plots.yaml"
    params:
        region_start = 6142346,
        region_end = 6177987
    shell:
        """
        python workflow/scripts/plots/heterozygosity_raw_plot.py \
            --input {input.tab} \
            --out_png {output.png} \
            --out_pdf {output.pdf} \
            --out_svg {output.svg} \
            --region_start {params.region_start} \
            --region_end {params.region_end} \
            > {log.out} 2> {log.err}
        """

################################################
## Rule: combined_heterozygosity_gene_plot
## Description: Combine smoothed heterozygosity plot with gene annotation track into one figure.
################################################

rule combined_heterozygosity_gene_plot:
    """
    Combine smoothed heterozygosity plot with gene annotation track into one figure.
    """
    input:
        vcf_tab = "tmp/amphioxus/amphioxus_chr4.tab",
        gff = "data/annotation/genomic.gff",
        top_snps = "results/snp/top5_snps_filtered.tsv"
    output:
        png = "results/plots/combined_view.png",
        pdf = "results/plots/combined_view.pdf",
        svg = "results/plots/combined_view.svg"
    log:
        out = "logs/plots/combined_plot.out",
        err = "logs/plots/combined_plot.err"
    conda:
        "../envs/plots.yaml"
    params:
        seqid = "OV696689.1",
        region_start = 6142346,
        region_end = 6177987
    shell:
        """
        python workflow/scripts/plots/combined_heterozygosity_gene_plot.py \
            --vcf_tab {input.vcf_tab} \
            --gff {input.gff} \
            --top_snps {input.top_snps} \
            --seqid {params.seqid} \
            --region_start {params.region_start} \
            --region_end {params.region_end} \
            --out_png {output.png} \
            --out_pdf {output.pdf} \
            --out_svg {output.svg} \
            > {log.out} 2> {log.err}
        """

################################################
## Rule: manhattan_plot
## Description: Create a Manhattan plot for LD clusters using Sex_g metric.
################################################

rule manhattan_sexg_plot:
    """
    Create a Manhattan plot for LD clusters using Sex_g metric.
    """
    input:
        candidates = "tmp/amphioxus/LD8.5cl20/candidates.csv"
    output:
        png = "results/plots/manhattan_sexg.png",
        pdf = "results/plots/manhattan_sexg.pdf",
        svg = "results/plots/manhattan_sexg.svg"
    log:
        out = "logs/plots/manhattan_sexg_plot.out",
        err = "logs/plots/manhattan_sexg_plot.err"
    conda:
        "../envs/plots.yaml"
    shell:
        """
        python workflow/scripts/plots/manhattan_sexg_plot.py \
            --input {input.candidates} \
            --out_png {output.png} \
            --out_pdf {output.pdf} \
            --out_svg {output.svg} \
            > {log.out} 2> {log.err}
        """

################################################
## Rule: manhattan_gc_adj_plot
## Description: Create a Manhattan plot for LD clusters using p_gc_adj metrics.
################################################

rule manhattan_gc_adj_plot:
    input:
        csv = "tmp/amphioxus/LD8.5cl20/candidates.csv"
    output:
        png = "results/plots/manhattan_gc_adj.png",
        pdf = "results/plots/manhattan_gc_adj.pdf",
        svg = "results/plots/manhattan_gc_adj.svg"
    log:
        out = "logs/plots/manhattan_gc_adj_plot.out",
        err = "logs/plots/manhattan_gc_adj_plot.err"
    conda:
        "../envs/plots.yaml"
    shell:
        """
        python workflow/scripts/plots/manhattan_gc_adj_plot.py \
            --input {input.csv} \
            --out_png {output.png} \
            --out_pdf {output.pdf} \
            --out_svg {output.svg} \
            > {log.out} 2> {log.err}
        """

rule manhattan_snp_plot:
    """
    Plot a Manhattan plot of sex-specific SNP association using genotype table.
    """
    input:
        tab = "tmp/amphioxus/amphioxus_all.tab"
    output:
        png = "results/plots/manhattan_snp.png",
        pdf = "results/plots/manhattan_snp.pdf",
        svg = "results/plots/manhattan_snp.svg",
        top_snps = "results/plots/top5_snps.json"
    log:
        out = "logs/plots/manhattan_snp.out",
        err = "logs/plots/manhattan_snp.err"
    conda:
        "../envs/plots.yaml"
    resources:
        mem_mb = 16000,
        cpus_per_task = 1,
        threads = 1,
        runtime = "2h"
    shell:
        """
        python workflow/scripts/plots/manhattan_snp_plot.py \
            --input {input.tab} \
            --out_png {output.png} \
            --out_pdf {output.pdf} \
            --out_svg {output.svg} \
            --out_top_snps {output.top_snps} \
            > {log.out} 2> {log.err}
        """


rule figure_2:
    """
    Combine smoothed heterozygosity plot with gene annotation track into one figure.
    """
    input:
        vcf_tab = "tmp/amphioxus/amphioxus_chr4.tab",
        gff = "data/annotation/genomic.gff"
    output:
        png = "results/plots/figure_2.png",
        pdf = "results/plots/figure_2.pdf",
        svg = "results/plots/figure_2.svg"
    log:
        out = "logs/plots/figure_2.out",
        err = "logs/plots/figure_2.err"
    conda:
        "../envs/plots.yaml"
    params:
        seqid = "OV696689.1",
        region_start = 6142346,
        region_end = 6177987
    shell:
        """
        python workflow/scripts/plots/figure_2.py \
            --vcf_tab {input.vcf_tab} \
            --gff {input.gff} \
            --seqid {params.seqid} \
            --region_start {params.region_start} \
            --region_end {params.region_end} \
            --out_png {output.png} \
            --out_pdf {output.pdf} \
            --out_svg {output.svg} \
            > {log.out} 2> {log.err}
        """

rule generate_combined_legend:
    output:
        png = "results/plots/legend_only.png",
        pdf = "results/plots/legend_only.pdf",
        svg = "results/plots/legend_only.svg"
    conda:
        "../envs/plots.yaml"
    shell:
        """
        python workflow/scripts/plots/generate_combined_legend.py \
            --out_png {output.png} \
            --out_pdf {output.pdf} \
            --out_svg {output.svg}
        """

rule heatmap_plot_flt1:
    input:
        "data/bgee/bgee_heatmap_flt1.tsv"
    output:
        png = "results/plots/heatmap_plot_flt1.png",
        pdf = "results/plots/heatmap_plot_flt1.pdf",
        svg = "results/plots/heatmap_plot_flt1.svg",
    log:
        out = "logs/plots/heatmap_plot_flt1.out",
        err = "logs/plots/heatmap_plot_flt1.err"
    conda:
        "../envs/plots.yaml"
    resources:
        mem_mb = 2000,
        cpus_per_task = 1,
        threads = 1,
        runtime = "2m"
    shell:
        """
        python workflow/scripts/plots/heatmap_plot.py \
            --input {input} \
            --out_png {output.png} \
            --out_pdf {output.pdf} \
            --out_svg {output.svg} \
            --gene FLT1 \
            > {log.out} 2> {log.err}
        """

rule heatmap_plot_hao1:
    input:
        "data/bgee/bgee_heatmap_hao1.tsv"
    output:
        png = "results/plots/heatmap_plot_hao1.png",
        pdf = "results/plots/heatmap_plot_hao1.pdf",
        svg = "results/plots/heatmap_plot_hao1.svg",
    log:
        out = "logs/plots/heatmap_plot_hao1.out",
        err = "logs/plots/heatmap_plot_hao1.err"
    conda:
        "../envs/plots.yaml"
    resources:
        mem_mb = 2000,
        cpus_per_task = 1,
        threads = 1,
        runtime = "2m"
    shell:
        """
        python workflow/scripts/plots/heatmap_plot.py \
            --input {input} \
            --out_png {output.png} \
            --out_pdf {output.pdf} \
            --out_svg {output.svg} \
            --gene HAO1 \
            > {log.out} 2> {log.err}
        """

rule heatmap_plot_dual:
    input:
        hao1 = "data/bgee/bgee_heatmap_hao1.tsv",
        flt1 = "data/bgee/bgee_heatmap_flt1.tsv"
    output:
        png = "results/plots/heatmap_plot_hao1_flt1.png",
        pdf = "results/plots/heatmap_plot_hao1_flt1.pdf",
        svg = "results/plots/heatmap_plot_hao1_flt1.svg",
    log:
        out = "logs/plots/heatmap_plot_hao1_flt1.out",
        err = "logs/plots/heatmap_plot_hao1_flt1.err"
    conda:
        "../envs/plots.yaml"
    resources:
        mem_mb = 2000,
        cpus_per_task = 1,
        threads = 1,
        runtime = "2m"
    shell:
        """
        python workflow/scripts/plots/heatmap_dual_plot.py \
            --input1 {input.hao1} --gene1 HAO1 \
            --input2 {input.flt1} --gene2 FLT1 \
            --out_png {output.png} \
            --out_pdf {output.pdf} \
            --out_svg {output.svg} \
            > {log.out} 2> {log.err}
        """

rule het_pc1_plot:
    input:
        rds = "tmp/amphioxus/LD8.5cl20/data_all.rds"
    output:
        png = "results/plots/het_pc1_plot.png",
        pdf = "results/plots/het_pc1_plot.pdf",
        svg = "results/plots/het_pc1_plot.svg"
    log:
        out = "logs/plots/het_pc1_plot.out",
        err = "logs/plots/het_pc1_plot.err"
    conda:
        "../envs/plots.yaml"
    resources:
        mem_mb = 2000,
        cpus_per_task = 1,
        threads = 1,
        runtime = "2m"
    shell:
        """
        Rscript workflow/scripts/plots/het_pc1_plot.R \
            --input {input.rds} \
            --out_png {output.png} \
            --out_pdf {output.pdf} \
            --out_svg {output.svg} \
            > {log.out} 2> {log.err}
        """

rule plot_normalized_ld_clusters:
    input:
        table="results/misc/normalized_ld_clusters.tsv",
        snps="results/misc/snp_counts/summary_snp_counts.tsv"
    output:
        png="results/plots/ld_clusters_per_mb.png",
        pdf="results/plots/ld_clusters_per_mb.pdf",
        svg="results/plots/ld_clusters_per_mb.svg"
    log:
        err="logs/plots/ld_clusters_per_mb.err"
    conda:
        "../envs/plots.yaml"
    shell:
        """
        python workflow/scripts/plots/normalized_ld_clusters_plot.py \
            --input {input.table} \
            --snp_table {input.snps} \
            --output results/plots/ld_clusters_per_mb \
            2> {log.err}
        """
