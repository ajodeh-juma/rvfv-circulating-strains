![rvfvphylo](docs/images/rvfvphylo_logo.png)

[![DOI](https://zenodo.org/badge/451463165.svg)](https://zenodo.org/badge/latestdoi/451463165)
[![Nextflow](https://img.shields.io/badge/nextflow-DSL2-brightgreen.svg)](https://www.nextflow.io/)
[![run with conda](https://img.shields.io/badge/run%20with-conda-3EB049.svg?logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1D355C.svg)](https://sylabs.io/singularity/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/ajodeh-juma/rvfv-circulating-strains/blob/master/LICENSE)
[![Twitter Follow](https://img.shields.io/twitter/follow/john_juma.svg?style=social)](https://twitter.com/john_juma)

## Table of contents

- [Introduction](#introduction)
- [Pipeline summary](#pipeline-summary)
- [Requirements](#requirements)
- [Installation](#installation)
  - [Option A — HPC with Singularity/Apptainer (recommended)](#option-a--hpc-with-singularityapptainer-recommended)
  - [Option B — Local installation with Conda](#option-b--local-installation-with-conda)
- [Running on an HPC cluster (SLURM)](#running-on-an-hpc-cluster-slurm)
- [Quick start](#quick-start)
- [Usage](#usage)
- [Parameters](#parameters)
- [Output](#output)
- [Analysis examples](#analysis-examples)
- [Credits](#credits)
- [Contributions and support](#contributions-and-support)
- [Citations](#citations)

## Introduction

**rvfvphylo** is a [Nextflow](https://www.nextflow.io) (DSL2) bioinformatics pipeline for the
comparative genetic and evolutionary characterization of circulating and vaccine strains of
Rift Valley fever virus (RVFV), across its three genomic segments (S, M and L).

Every process is containerized, so the underlying tools (MAFFT, IQ-TREE, HyPhy, TreeTime,
BCFtools, SnpSift, snipit, seqkit, 3SEQ, TrimAl, snp-sites and a set of R/Python helper
scripts) never need to be installed by hand — Nextflow pulls a per-process container image
automatically at run time when using the `docker`, `singularity` or (on most HPC systems)
`apptainer` engine. A `conda` environment is provided as a fallback for systems where
containers are not available.

## Pipeline summary

The pipeline is organized into three independently runnable subworkflows, selected with
`--subworkflow`:

1. **`rvfvcirculatingstrains`** — comparative genetic and evolutionary analysis of circulating
   sequences for a given segment against a reference vaccine strain (Smithburn, MP-12 or
   Clone 13). Produces alignments, a maximum-likelihood phylogeny (IQ-TREE), selection
   analysis (HyPhy FEL/FUBAR/SLAC), strain typing/SNP calls, and a time-calibrated
   phylogeography (TreeTime).
2. **`rvfvmutationalprofiling`** — mutational profiling of circulating sequences against the
   ZH548 reference strain, including host-associated sequence-context (APOBEC-like editing
   signature) analysis.
3. **`rvfvphylocontinuous`** — continuous phylogeographic analysis (geocoding, divergence
   dating) of a segment's sequence set.

## Requirements

To **run** the pipeline you need, at minimum:

- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html) `>=19.10.0` (developed and
  tested with `23.04`)
- One of: [Singularity](https://docs.sylabs.io/guides/latest/user-guide/)/[Apptainer](https://apptainer.org/)
  (recommended, and typically what's available on HPC clusters), Docker, or Conda/Mamba
- Java 11+ (required by Nextflow)

You do **not** need `conda`, `bioconda`, or `singularity`/`apptainer` themselves to be
buildable or installable on your laptop — they only need to be usable on whichever machine
actually executes the pipeline (see below).

## Installation

```bash
git clone https://github.com/ajodeh-juma/rvfv-circulating-strains.git
cd rvfv-circulating-strains
```

How you go from here depends on where the pipeline will actually run.

### Option A — HPC with Singularity/Apptainer (recommended)

This is the supported way to run the pipeline when your local machine (e.g. macOS, or any
system without root/Singularity support) cannot install Conda environments or build/run
Singularity containers. All of the compute happens on the HPC cluster; your laptop is only
used to edit files and submit the job.

1. Copy (or `git clone` directly on) the cluster:

   ```bash
   rsync -avP --exclude=".git" --exclude="work" --exclude="output-dir" \
     ./rvfv-circulating-strains <user>@<hpc-login-node>:~/projects/pipelines/
   ```

2. Log in to the HPC and load (or install) Nextflow and Singularity/Apptainer. Cluster module
   names vary by site — check with `module avail nextflow singularity apptainer`, e.g.:

   ```bash
   module load nextflow
   module load singularity   # or: module load apptainer
   ```

   If no `nextflow` module exists, install it into a lightweight conda env (only Nextflow
   itself, none of the pipeline's tools):

   ```bash
   conda create -n nextflow -c bioconda nextflow=23.04
   conda activate nextflow
   ```

3. Run with `-profile singularity`. On first run, Nextflow will automatically pull the
   required per-process Singularity images from
   [Galaxy Depot](https://depot.galaxyproject.org/singularity/) /
   [Biocontainers](https://quay.io/organization/biocontainers) and cache them under
   `work/singularity` — no manual container build step is needed:

   ```bash
   nextflow run main.nf -profile singularity --help
   ```

See [Running on an HPC cluster (SLURM)](#running-on-an-hpc-cluster-slurm) for a full example
submission workflow.

### Option B — Local installation with Conda

Use this only on a system where Conda/Bioconda packages install and run natively (Linux is
best supported; some tools in `environment.yml` are not available for macOS/arm64).

```bash
conda env create -n rvfv-phylo -f environment.yml
conda activate rvfv-phylo
nextflow run main.nf -profile conda --help
```

## Running on an HPC cluster (SLURM)

Nextflow submits and supervises each process as its own cluster job, so the `nextflow run`
command itself must keep running for the whole pipeline duration. On a SLURM cluster there
are two common ways to do that:

**1. Run Nextflow from a persistent session on the login node**, letting it submit
`sbatch` jobs for each process:

```bash
# on the login node
screen -S rvfvphylo   # or: tmux new -s rvfvphylo

module load nextflow singularity
export NXF_OPTS='-Xms1g -Xmx4g'

nextflow run main.nf \
    -profile singularity \
    --subworkflow rvfvcirculatingstrains \
    --fasta data/RVFV-M.fasta \
    --metadata data/RVFV-M.csv \
    --segment M \
    --prefix riftM \
    --start 20 --end 3611 \
    --vaccine_reference DQ380208 \
    --outliers assets/RVFV-M-outliers-circulating.txt \
    --outdir output-dir/riftM \
    -work-dir work/riftM \
    -resume

# detach with Ctrl-A D (screen) / Ctrl-B D (tmux); check back later with:
#   screen -r rvfvphylo   or   tmux attach -t rvfvphylo
```

This requires an `executor { name = 'slurm' }` (or equivalent) entry so that individual
processes are submitted to the scheduler instead of running on the login node — add a
site-specific config, e.g. `conf/hpc.config`:

```groovy
process {
  executor = 'slurm'
  queue    = 'batch'          // replace with your cluster's partition/queue name
}
executor {
  queueSize = 20
}
```

and include it on the command line with `-c conf/hpc.config` (or add it to `nextflow.config`
under its own profile) alongside `-profile singularity`.

**2. Submit the whole Nextflow head process as a single batch job** (simplest, but ties up
one job allocation for the full pipeline runtime):

```bash
#!/bin/bash
#SBATCH --job-name=rvfvphylo
#SBATCH --output=rvfvphylo.%j.log
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=48:00:00

module load nextflow singularity

nextflow run main.nf -profile singularity -c conf/hpc.config \
    --subworkflow rvfvcirculatingstrains \
    --fasta data/RVFV-M.fasta \
    --metadata data/RVFV-M.csv \
    --segment M --prefix riftM --start 20 --end 3611 \
    --vaccine_reference DQ380208 \
    --outliers assets/RVFV-M-outliers-circulating.txt \
    --outdir output-dir/riftM -work-dir work/riftM -resume
```

submitted with `sbatch run_rvfvphylo.sbatch`. In both cases, `-resume` lets you re-launch the
same command after a failure or a walltime cut-off without recomputing already-finished
steps.

## Quick start

Show all pipeline options:

```bash
nextflow run main.nf --help
```

Show all options, including advanced/hidden ones:

```bash
nextflow run main.nf --help --show_hidden_params
```

## Usage

1. **`rvfvcirculatingstrains`** — genetic and phylogenetic analysis of the M segment against
   the MP-12 vaccine strain:

   ```bash
   BASEDIR="${HOME}/projects/vaccine-and-circulating-strains-analysis/segments/M/complete/global"

   nextflow run main.nf \
       -profile singularity \
       --subworkflow rvfvcirculatingstrains \
       --fasta ${BASEDIR}/merged-sequences/RVFV-M.fasta \
       --metadata ${BASEDIR}/merged-sequences/RVFV-M.csv \
       --segment M \
       --prefix riftM \
       --start 20 \
       --end 3611 \
       --remove_duplicates false \
       --lineages ${BASEDIR}/assignment/output-dir/report/lineages.csv \
       --vaccine_reference DQ380208 \
       --outliers ./assets/RVFV-M-outliers-circulating.txt \
       --recombinants ./assets/rvfv-M-potential-recombinants.txt \
       --outdir ./output-dir/riftM \
       -work-dir ./work/riftM \
       -resume
   ```

2. **`rvfvmutationalprofiling`** — should follow step 1. Profiles mutations against the
   ZH548 reference strain:

   ```bash
   nextflow run main.nf \
       -profile singularity \
       --subworkflow rvfvmutationalprofiling \
       --fasta ${BASEDIR}/merged-sequences/rvfv-M.fasta \
       --metadata ${BASEDIR}/merged-sequences/rvfv-M.csv \
       --segment M \
       --prefix riftM \
       --start 20 \
       --end 3611 \
       --remove_duplicates false \
       --lineages ${BASEDIR}/assignment/output-dir/report/lineages.csv \
       --vaccine_reference NC_014396 \
       --outliers ./assets/RVFV-M-outliers-circulating.txt \
       --recombinants ./assets/rvfv-M-potential-recombinants.txt \
       --outdir ./output-dir/mutational-profiling/riftM \
       -work-dir ./work/mutational-profiling/riftM \
       -resume
   ```

   Then use `bin/getSequenceContexts.py` to extract sequence contexts from the alignments,
   and `scripts/plots.R` to render the summary figures.

Swap `-profile singularity` for `-profile conda` if running under Option B, or
`-profile docker` if Docker is available instead.

## Parameters

```text
Input/output options
  --subworkflow        [string]  'rvfvcirculatingstrains', 'rvfvmutationalprofiling' or 'rvfvphylocontinuous'
  --fasta              [string]  Input FASTA file containing the sequences
  --metadata           [string]  Input CSV metadata file (requires 'sample_name' and 'Ct' columns)
  --recombinants       [string]  Input text file of recombinant sequence accessions (column 'recombinants')
  --segment            [string]  Genomic segment: 'S-NSS', 'S-NP', 'M' or 'L'
  --outdir             [string]  Output directory [default: ./results]
  --email              [string]  Email address for completion summary

Alignment masking options
  --start              [integer] Start position to trim the alignment (0-based)
  --end                [integer] End position to trim the alignment (0-based)

Alignment filtering options
  --outliers           [string]  Text file of outlier sequence identifiers (reformatted headers, as
                                  generated by the REFORMAT_HEADERS process)

Strain characterization options
  --lineages           [string]  CSV lineage-assignment file (as generated by the rvfvtyper pipeline)
  --vaccine_reference  [string]  Vaccine-strain reference accession (must be present in your sequence set)
  --grouping_column    [string]  Column used to group sequences [default: strain_type]
  --min_freq           [number]  Minimum frequency to support a SNP call [default: 0.2]
  --max_freq           [number]  Maximum frequency to support a SNP call [default: 0.8]
  --seq_type           [string]  'dna' or 'protein' [default: protein]
  --snp_type           [string]  'singleton', 'multiple' or 'conserved' [default: singleton]
  --group_per_lineage  [string]  Group the stats output per lineage [default: false]

Dataset filtering options
  --filter_columns     [string]  Space-separated metadata columns required to be populated
                                  (e.g. 'country location host date') [default: location]

Process skipping options
  --skip_modeltesting  [boolean] Skip the modeltest-ng model-selection step

Execution profiles (-profile)
  singularity                    Run every process in its own Singularity/Apptainer container (recommended for HPC)
  docker                         Run every process in its own Docker container
  conda                          Create/use a single Conda environment from environment.yml
  test                           Minimal test configuration
```

Run `nextflow run main.nf --help --show_hidden_params` for the complete, always up-to-date
list (sourced from [`nextflow_schema.json`](nextflow_schema.json)).

## Output

For **`rvfvcirculatingstrains`**, the typical outputs are laid out as:

```text
output-dir/riftM/
├── alignment
│   ├── riftM.log
│   ├── riftM_align.fasta
│   ├── riftM_dedup.fasta
│   ├── riftM_dedup.txt
│   ├── riftM_duplicated_taxa.txt
│   ├── riftM_filtered.fasta
│   ├── riftM_filtered.txt
│   └── riftM_masked.fasta
├── hyphy
│   ├── riftM.FEL.json
│   ├── riftM.FEL.log
│   ├── riftM.FUBAR.cache
│   ├── riftM.FUBAR.json
│   ├── riftM.FUBAR.log
│   ├── riftM.SLAC.json
│   └── riftM.SLAC.log
├── iqtree
│   ├── riftM.bionj
│   ├── riftM.iqtree
│   └── riftM.treefile
├── models
│   ├── riftM.model.ckp
│   ├── riftM.model.log
│   ├── riftM.model.out
│   ├── riftM.model.topos
│   └── riftM.model.tree
├── phylogeo
│   ├── riftM_dates.csv
│   ├── riftM_geocoded.pdf
│   ├── riftM_geocoded.txt
│   └── riftM_geolocations.csv
├── pipeline_info
│   ├── execution_report.html
│   ├── execution_timeline.html
│   ├── execution_trace.txt
│   └── pipeline_dag.svg
├── reformatted-sequences
│   ├── riftM.fasta
│   └── riftM.txt
├── sequences
│   ├── riftM.csv
│   └── riftM.fasta
├── strain-types
│   ├── riftM.DQ380208.all.csv
│   ├── riftM.DQ380208.per.lineage.csv
│   ├── riftM.DQ380208.snps.csv
│   ├── riftM.DQ380208.snps.pdf
│   ├── riftM.DQ380208.txt
│   ├── riftM.fasta
│   ├── riftM.strain_type.DQ380208.amino.acid.fasta
│   ├── riftM.strain_type.DQ380208.bcftools.stats.tstv.txt
│   ├── riftM.strain_type.DQ380208.mutations.per.strain.singleton.csv
│   ├── riftM.strain_type.DQ380208.parsed.bcftools.stats.csv
│   ├── riftM.strain_type.DQ380208.singleton.txt
│   ├── riftM.strain_type.DQ380208.snpSift.tstv.txt
│   ├── riftM.strain_type.DQ380208.vcf
│   ├── riftM.traits.txt
│   ├── strain_type.DQ380208.sorted.alignment.all.labels.csv
│   └── strain_type.DQ380208.sorted.alignment.fasta
└── treetime
    ├── divergence_tree.nexus
    ├── riftM.log
    ├── riftM_ancestral_sequences.fasta
    ├── riftM_timetree.nexus
    └── trace_run.log
```

For **`rvfvmutationalprofiling`**, the typical outputs are laid out as:

```text
output-dir/mutational-profiling/riftM/
├── alignment
│   ├── riftM.log
│   ├── riftM.reverse.complement.fasta
│   ├── riftM_align.fasta
│   ├── riftM_dedup.fasta
│   ├── riftM_dedup.txt
│   ├── riftM_duplicated_taxa.txt
│   ├── riftM_filtered.fasta
│   ├── riftM_filtered.txt
│   └── riftM_masked.fasta
├── pipeline_info
│   ├── execution_report.html
│   ├── execution_timeline.html
│   ├── execution_trace.txt
│   └── pipeline_dag.svg
├── reformatted-sequences
│   ├── riftM.fasta
│   └── riftM.txt
├── sequence-contexts
│   ├── A3A_A3B.CT.riftM.context.txt
│   ├── A3A_A3B.GA.riftM.context.txt
│   ├── A3A_A3B.riftM.edited.sites.txt
│   ├── ... (per APOBEC-family editor, per host species: antelope, bat, buffalo, cow, human,
│   │        mosquito, sheep)
│   ├── riftM-CT-sequence-contexts-with-metadata.csv
│   ├── riftM-GA-sequence-contexts-with-metadata.csv
│   └── riftM-sequence-contexts-with-metadata.csv
├── sequences
│   ├── riftM.csv
│   └── riftM.fasta
├── strain-types
│   ├── riftM.NC_014396.txt
│   ├── riftM.fasta
│   ├── riftM.strain_type.NC_014396.mutations.per.strain.singleton.csv
│   ├── riftM.strain_type.NC_014396.singleton.txt
│   ├── riftM.strain_type.NC_014396.vcf
│   ├── riftM.traits.txt
│   ├── strain_type.NC_014396.sorted.alignment.all.labels.csv
│   └── strain_type.NC_014396.sorted.alignment.fasta
└── strain-types-rev
    └── (mirrors strain-types/, computed on the reverse-complement strand)
```

Use `bin/getSequenceContexts.py` to extract sequence contexts from the alignments, and
`scripts/plots.R` to generate the summary figures from them.

## Analysis examples

[`analysis/vaccine-and-circulating-strains/`](analysis/vaccine-and-circulating-strains/README.md)
is a worked, dataset-specific example that post-processes per-segment
`rvfvcirculatingstrains` runs (root-to-tip plots, strain-type/lineage summaries, vaccine
substitution tables) into manuscript-ready figures and tables.

## Credits

rvfvphylo was originally written by [ajodeh-juma](https://github.com/ajodeh-juma).

## Contributions and support

If you would like to contribute to this pipeline, please see the
[contributing guidelines](.github/CONTRIBUTING.md). Issues and feature requests can be filed
on the [issue tracker](https://github.com/ajodeh-juma/rvfv-circulating-strains/issues).

## Citations

<!-- If you use rvfvphylo for your analysis, please cite it using the following doi:
[10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

This pipeline was built from the [nf-core](https://nf-co.re) pipeline template. If you use
tooling from that template, please cite:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas
> Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi:
> [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
> ReadCube: [Full Access Link](https://rdcu.be/b1GjZ)
