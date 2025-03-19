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
