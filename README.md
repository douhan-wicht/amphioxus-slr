# The Foo Project

This project serves as a sandbox for familiarizing myself with the directory structure I plan to use in my master project, which will extensively leverage Snakemake.

## Basic Environment

- **conda**: *Insert version here*
- **mamba**: *Insert version here*
- **snakemake**: 8.29.0

---

## Before Starting

Activate the Snakemake environment:

```sh
micromamba activate snakemake-env
```

---

## Snakemake Usage

### Dry Run

For an initial dry run, execute:

```sh
snakemake --cores 1 -p --use-conda -n
```

### Cluster Run

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

---

## Folder Organization

- `config`  → Configuration files (modify for different genomes or datasets).
- `rules`  → Snakemake rules calling scripts to generate results.
- `scripts`  → Scripts for raw data analysis and result generation.
- `envs`  → Conda environments corresponding to specific rule groups.
- `metadata`  → Sample metadata.
- `results`  → Intermediate results & plots (organized by rule groups).
- `data`  → Raw data (FASTQ files, genome files, etc.).
- `logs`  → Logs for debugging and tracking runs.
- `benchmarks`  → Performance benchmarks of workflows.

---

## Workflow

This section should describe the workflow execution order, dependencies, and outputs. Include details on:

- Which rules are executed first.
- What files are produced at each step.
- The reasoning behind the structure.
