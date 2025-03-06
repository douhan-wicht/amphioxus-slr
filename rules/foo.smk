###########################################################################

#  ███████╗ ██████╗  ██████╗ 
#  ██╔════╝██╔═══██╗██╔═══██╗
#  ██████╗ ██║   ██║██║   ██║
#  ██╔══╝  ██║   ██║██║   ██║
#  ██║     ╚██████╔╝╚██████╔╝
#  ╚═╝      ╚═════╝  ╚═════╝ 
#
# This file contains rules for processing and analyzing data related to the "foo" module.
# It includes rules for generating specific outputs, such as plots or processed data files,
# and is part of a larger Snakemake workflow for reproducing the paper's results.

###########################################################################

################################################
## Rule: foo_1
## Description: This rule creates foo 1 ...
################################################

rule foo_1:
    '''
    This rule creates foo 1 ...
    '''
    input:
        Metadata = "metadata/foo_metadata.txt",
    output:
        PDF = "results/foo/foo_1.pdf",
    log:
        err = "logs/foo/foo_1.err",
        out = "logs/foo/foo_1.out"
    benchmark:
        "benchmarks/foo/foo_1.txt"
    conda:
        '../envs/foo.yaml'
    params:
        time = '1:00:00',
        name = "pFig1",
        threads = 1,
        mem = 2000,
    shell:
        """
        scripts/foo/foo_1.R > {log.out} 2> {log.err}
        """
