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
        directory("data/DNAseqVCF")
    log:
        err = "logs/setup/import_data.err",
        out = "logs/setup/import_data.out"
    conda:
        "../envs/setup.yaml"
    params:
        name = "import_data",
        time = "00:30:00"
    resources:
        mem = 2000,
        threads = 1
    shell:
        """
        mkdir -p {output}
        rsync -av --progress /nas/FAC/FBM/DEE/mrobinso/default/D2c/mbrasovi/Banyuls_Roscoff/DNAseqVCF/ {output} > {log.out} 2> {log.err}
        """
