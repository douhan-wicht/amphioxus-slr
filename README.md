# The Amphioxus SLR Project

![Project Logo](ressources/amphioxus_logo.png)

## Overview

The **Amphioxus SLR Project** aims to apply a novel tool, [SLRfinder](https://github.com/xuelingyi/SLRfinder), to identify sex-linked regions in *Amphioxus*. By leveraging this tool, we seek to gain insights into the mechanisms of sex determination in *Amphioxus* and contribute to a broader understanding of the diversification of sex-determination mechanisms in chordates.

## Environment Requirements

Ensure you have the following dependencies installed:

- **micromamba**: *Insert version here*
- **conda**: *Insert version here*
- **snakemake**: 8.29.0

## Setup

Before running the workflow, activate the Snakemake environment:

```sh
micromamba activate snakemake-env
```

## Running Snakemake

### Dry Run (Validation)

To perform a dry run and validate the workflow without executing commands:

```sh
snakemake --cores 1 -p --use-conda -n
```

### Cluster Execution

For execution on a computing cluster:

```sh
snakemake -p -j 30 \
    --snakefile workflow/snakefile \
    --directory . \
    --cluster "sbatch -J {params.name} -N 1 \
    -o logs/.slurm/%x.out -e logs/.slurm/%x.err \
    --cpus-per-task={params.threads} \
    --mem={params.mem} -t {params.time}" \
    --use-conda
```

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

## Workflow Overview

This section provides details on the workflow execution, dependencies, and outputs:

- **Execution Order:** Describes the sequence in which rules are executed.
- **Generated Files:** Specifies the outputs at each stage.
- **Workflow Rationale:** Explains the design choices and dependencies.

_(Consider expanding this section with a visual workflow diagram or example outputs.)_

---

For any questions or contributions, please refer to the project's documentation or contact the maintainers.
