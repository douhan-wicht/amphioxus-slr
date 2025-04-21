rule annotate_snps:
    """
    Converts SNPs from JSON to BED format, then uses bedtools to annotate them using a GTF file.
    """
    input:
        json = "results/plots/top5_snps.json",
        annotation = "data/annotation/genomic.gtf"
    output:
        annotated = "results/snp/top5_snps_annotated.tsv"
    log:
        err = "logs/snp/annotate_snps.err",
        out = "logs/snp/annotate_snps.out"
    conda:
        "../envs/snp.yaml"
    resources:
        mem_mb = 2000,
        cpus_per_task = 1,
        threads = 1,
        runtime = "10m"
    shell:
        """
        TMP_BED=$(mktemp)

        python workflow/scripts/snp/json_to_bed.py \
            --json {input.json} \
            --bed $TMP_BED \
            >> {log.out} 2>> {log.err}

        bedtools intersect \
            -a $TMP_BED \
            -b {input.annotation} \
            -wa -wb \
            > {output.annotated} \
            2>> {log.err}
        """

rule extract_minimal_gff:
    input:
        gff="data/annotation/genomic.gff"
    output:
        mini_gff="data/annotation/FLT1_HAO1_subset.gff"
    params:
        loci=["BLAG_LOCUS17194", "BLAG_LOCUS17195"]
    conda:
        "../envs/snp.yaml"
    shell:
        """
        python workflow/scripts/snp/extract_locus_gff.py \
            --input {input.gff} \
            --output {output.mini_gff} \
            --loci {params.loci}
        """

rule fix_gff_for_snpeff:
    input:
        "data/annotation/FLT1_HAO1_subset.gff"
    output:
        "data/annotation/FLT1_HAO1_subset.snpeff.gff"
    conda:
        "../envs/snp.yaml"
    shell:
        """
        python workflow/scripts/snp/fix_gff_for_snpeff.py --input {input} --output {output}
        """


rule build_snpeff_db:
    """
    Creates a custom SnpEff DB for Branchiostoma using files in data/annotation/.
    This rule ensures that the required snpEff.config is created and properly configured.
    """
    input:
        fasta = "data/annotation/GCA_927797965.1_BraLan3_genomic.fna",
        gff = "data/annotation/FLT1_HAO1_subset.snpeff.gff" # Use mini or full gff
    output:
        db_built_flag = "data/annotation/snpeff/BranchiostomaLanceolatum/build.done"
    log:
        "logs/snp/build_snpeff.err"
    conda:
        "../envs/snp.yaml"
    shell:
        r"""
        set -euo pipefail

        echo "Creating SnpEff database folder..."
        mkdir -p data/annotation/snpeff/BranchiostomaLanceolatum

        echo "Linking genome FASTA..."
        ln -sf $(realpath {input.fasta}) data/annotation/snpeff/BranchiostomaLanceolatum/sequences.fa

        echo "Linking GFF..."
        ln -sf $(realpath {input.gff}) data/annotation/snpeff/BranchiostomaLanceolatum/genes.gff

        echo "Writing snpEff.config..."
        DATA_DIR=$(realpath data/annotation/snpeff)
        cat > data/annotation/snpeff/snpEff.config << EOF
data.dir = ${{DATA_DIR}}
BranchiostomaLanceolatum.genome : BranchiostomaLanceolatum
EOF

        echo "Building SnpEff database..."
        snpEff build -gff3 -noCheckCds -noCheckProtein -c data/annotation/snpeff/snpEff.config BranchiostomaLanceolatum

        echo "Done building SnpEff DB."
        touch {output.db_built_flag}
        """

rule annotate_snps_with_snpeff:
    """
    Annotate chr4 variants using the custom SnpEff DB.
    """
    input:
        vcf = "tmp/amphioxus/a15m75/amphioxus_chr4_a15m75.recode.vcf",
        db_flag = rules.build_snpeff_db.output.db_built_flag
    output:
        vcf = "results/snp/chr4_snps.ann.vcf"
    log:
        "logs/snp/annotate_snps.err"
    conda:
        "../envs/snp.yaml"
    shell:
        """
        snpEff -c data/annotation/snpeff/snpEff.config \
            -v BranchiostomaLanceolatum \
            {input.vcf} > {output.vcf} 2>> {log}
        """

rule summarize_chr4_slrs:
    input:
        vcf="results/snp/chr4_snps.ann.vcf"
    output:
        summary="results/snp/chr4_hao1_flt1_snp_summary.tsv"
    conda:
        "../envs/snp.yaml"
    shell:
        """
        python workflow/scripts/snp/parse_snpeff_chr4.py --input {input.vcf} --output {output.summary}
        """

rule filter_top_snps:
    input:
        json="results/plots/top5_snps.json",
        snps="tmp/amphioxus/amphioxus_chr4.tab"
    output:
        top_snps="results/snp/top5_snps_filtered.tsv"
    conda:
        "../envs/snp.yaml"
    shell:
        """
        python workflow/scripts/snp/filter_top_snps.py \
            --json {input.json} \
            --snps {input.snps} \
            --output {output.top_snps}
        """

rule classify_snps:
    input:
        snps="results/snp/top5_snps_filtered.tsv",
        ref="data/annotation/GCA_927797965.1_BraLan3_genomic.fna",
        gff="data/annotation/genomic.gff"
    output:
        effects="results/snp/top5_snp_effects.tsv"
    conda:
        "../envs/snp.yaml"
    shell:
        """
        python workflow/scripts/snp/classify_snps.py \
            --snps {input.snps} \
            --ref {input.ref} \
            --gff {input.gff} \
            --output {output.effects}
        """