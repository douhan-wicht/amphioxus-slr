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
    # conda:
    #     "/work/FAC/FBM/DEE/mrobinso/evolseq/dwicht1/envs/setup/setup.yaml"
    params:
        time = "1:00:00",
        name = "import_data",
        threads = 1,
        mem = "2G",
    shell:
        """
        mkdir -p {output}
        rsync -av --progress /nas/FAC/FBM/DEE/mrobinso/default/D2c/mbrasovi/Banyuls_Roscoff/DNAseqVCF/ {output} > {log.out} 2> {log.err}
        """
