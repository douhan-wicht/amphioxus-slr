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
        "results/plots/karyotype.png",
        "results/plots/karyotype.pdf"
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
            --out_png {output[0]} \
            --out_pdf {output[1]} \
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
        "results/plots/gene_region.png",
        "results/plots/gene_region.pdf"
    params:
        seqid = "OV696689.1",
        start = 6142346,
        end = 6164195
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
            --out_png {output[0]} \
            --out_pdf {output[1]} \
            > {log.out} 2> {log.err}
        """

################################################
## Rule: vcf_to_tab
## Description: Convert VCF to tabular format using GATK VariantsToTable.
################################################

# rule vcf_to_tab:
#     """
#     Convert VCF to tabular format using GATK VariantsToTable.
#     """
#     input:
#         vcf = "tmp/amphioxus/a15m75/amphioxus_chr4_a15m75.recode.vcf"
#     output:
#         tab = "tmp/amphioxus/amphioxus_chr4.tab"
#     log:
#         out = "logs/plots/vcf_to_tab.out",
#         err = "logs/plots/vcf_to_tab.err"
#     conda:
#         "../envs/plots.yaml"
#     shell:
#         """
#         gatk VariantsToTable \
#             -V {input.vcf} \
#             -O {output.tab} \
#             -F CHROM -F POS -F REF -F ALT \
#             -GF GT \
#             > {log.out} 2> {log.err}
#         """

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
        pdf = "results/plots/heterozygosity_plot.pdf"
    log:
        out = "logs/plots/heterozygosity_plot.out",
        err = "logs/plots/heterozygosity_plot.err"
    conda:
        "../envs/plots.yaml"
    params:
        region_start = 6142346,
        region_end = 6164195
    shell:
        """
        python workflow/scripts/plots/heterozygosity_plot.py \
            --input {input.tab} \
            --out_png {output.png} \
            --out_pdf {output.pdf} \
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
        pdf = "results/plots/heterozygosity_raw.pdf"
    log:
        out = "logs/plots/heterozygosity_raw_plot.out",
        err = "logs/plots/heterozygosity_raw_plot.err"
    conda:
        "../envs/plots.yaml"
    params:
        region_start = 6142346,
        region_end = 6164195
    shell:
        """
        python workflow/scripts/plots/heterozygosity_raw_plot.py \
            --input {input.tab} \
            --out_png {output.png} \
            --out_pdf {output.pdf} \
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
        gff = "data/annotation/genomic.gff"
    output:
        png = "results/plots/combined_view.png",
        pdf = "results/plots/combined_view.pdf"
    log:
        out = "logs/plots/combined_plot.out",
        err = "logs/plots/combined_plot.err"
    conda:
        "../envs/plots.yaml"
    params:
        seqid = "OV696689.1",
        region_start = 6142346,
        region_end = 6164195
    shell:
        """
        python workflow/scripts/plots/combined_heterozygosity_gene_plot.py \
            --vcf_tab {input.vcf_tab} \
            --gff {input.gff} \
            --seqid {params.seqid} \
            --region_start {params.region_start} \
            --region_end {params.region_end} \
            --out_png {output.png} \
            --out_pdf {output.pdf} \
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
        pdf = "results/plots/manhattan_sexg.pdf"
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
        pdf = "results/plots/manhattan_gc_adj.pdf"
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
        pdf = "results/plots/manhattan_snp.pdf"
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
            > {log.out} 2> {log.err}
        """
