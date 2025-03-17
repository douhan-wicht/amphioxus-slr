# The Amphioxus SLR Project

![Project Logo](ressources/amphioxus_logo.png)

## Overview

The **Amphioxus SLR Project** aims to apply a novel tool, [SLRfinder](https://github.com/xuelingyi/SLRfinder), to identify sex-linked regions in *Amphioxus*. By leveraging this tool, we seek to gain insights into the mechanisms of sex determination in *Amphioxus* and contribute to a broader understanding of the diversification of sex-determination mechanisms in chordates.

## Environment Requirements

Ensure you have the following dependencies installed:

- **conda**: 25.1.1
- **snakemake**: 8.30.0

## Setup

Before running the workflow, activate the Snakemake environment containing snakemake and conda:

```sh
conda activate snakemake-8.30.0
```

The workflow was setup according to the best practices described on the snakemake Github repository.

## Running Snakemake

### Dry Run (Validation)

To perform a dry run and validate the workflow without executing commands:

```sh
snakemake --cores 1 -p --use-conda -n
```

### Cluster Execution

The workflow can then be executed using:
```sh
snakemake -p -j 30 --profile curnagl --use-conda --executor slurm
```

STILL A LOT TO CONFIGURE FOR THE snakemake profile so that the slurm executor works properly. Go see documentation on https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html.


## Project Directory Structure

The project follows a structured organization to facilitate reproducibility and ease of use:

- **`config/`**  → Configuration files (modify for different genomes or datasets).
- **`rules/`**  → Snakemake rules defining the workflow.
- **`scripts/`**  → Scripts for raw data analysis and result generation.
- **`envs/`**  → Conda environments for specific workflow steps.
- **`metadata/`**  → Sample metadata.
- **`results/`**  → Intermediate and final results, including plots (organized by rule groups).
- **`data/`**  → Raw data files (e.g., FASTQ files, genome sequences).
- **`logs/`**  → Log files for debugging and tracking runs.
- **`benchmarks/`**  → Workflow performance benchmarks.

The project directory structure is as follows:

```plaintext
├── README.md
├── benchmarks
├── config
│   └── amphioxus-slr.yaml
├── data
│   └── DNAseqVCF
├── envs
│   ├── SLRfinder.yaml
│   └── setup.yaml
├── logs
│   └── setup
│       ├── import_data.err
│       └── import_data.out
├── metadata
│   ├── metadata.csv
│   └── reference.list
├── ressources
│   └── amphioxus_logo.png
├── results
│   └── SLRfinder
├── rules
│   ├── SLRfinder.smk
│   └── setup.smk
├── scripts
│   ├── SLRfinder
│   │   └── amphioxus
│   │       └── SLRfinder_functions.r
│   └── setup
└── snakefile
```

## Workflow Overview

This section provides details on the workflow execution, dependencies, and outputs:

- **Execution Order:** Describes the sequence in which rules are executed.
- **Generated Files:** Specifies the outputs at each stage.
- **Workflow Rationale:** Explains the design choices and dependencies.
