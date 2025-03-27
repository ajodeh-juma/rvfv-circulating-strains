<!-- # ![rvfvphylo](docs/images/rvfvphylo_logo.png) -->


- [Introduction](#introduction)
  - [Installation](#installation)
  - [Usage](#usage)
  - [Method details](#method-details)
  - [Output](#output)
  - [Pipeline summary](#pipeline-summary)
  - [Citations](#citations)
  


[![DOI](https://zenodo.org/badge/451463165.svg)](https://zenodo.org/badge/latestdoi/451463165)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A520.04.0-brightgreen.svg)](https://www.nextflow.io/)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](https://bioconda.github.io/)
[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23rvfvtyping-4A154B?logo=slack)](https://nfcore.slack.com/channels/rvfvtyping)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-brightgreen.svg)](https://github.com/ajodeh-juma/rvfvtyping/blob/master/LICENSE)
<!-- [![Docker](https://img.shields.io/docker/automated/nfcore/rvfvtyping.svg)](https://hub.docker.com/r/nfcore/rvfvtyping) -->
<!-- [![GitHub Actions CI Status](https://github.com/nf-core/rvfvtyping/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/rvfvtyping/actions) -->
<!-- [![GitHub Actions Linting Status](https://github.com/nf-core/rvfvtyping/workflows/nf-core%20linting/badge.svg)](https://github.com/nf-core/rvfvtyping/actions) -->
<!-- ![](https://img.shields.io/badge/uses-docker-blue.svg) -->

[![Twitter
Follow](https://img.shields.io/twitter/follow/john_juma.svg?style=social)](https://twitter.com/john_juma)


## Introduction

<!-- TODO nf-core: Add a brief overview of what the pipeline does and how it works -->

**rvfvphylo** is a bioinformatics pipeline for characterizing circulating and vaccine strains
of the Rift Valley fever virus.
The pipeline has 3 main subworkflows:

1. **rvfvcirculatingstrains**: This subworkflow performs a comparative genetic and evolutionary
   analysis for the 3 segments of RVFV against the commonly used vaccine strains
   (Smithburn, MP-12 and Clone-13)
2. **rvfvmutationalprofiling**: This subworkflow performs mutational profiling for
   the 3 segments using the ZH548 strain.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool
to run tasks across multiple compute infrastructures in a very portable manner.
It comes with docker containers making installation trivial and results highly
reproducible.

## Installation

**rvfvphylo** runs on UNIX/LINUX systems. You will install Miniconda3 from [here](https://docs.conda.io/en/latest/miniconda.html). Once Miniconda3 has been installed, proceed with pipeline installation

```
git clone https://github.com/ajodeh-juma/rvfv-circulating-strains.git
cd rvfv-circulating-strains
conda env create -n rvfv-phylo -f environment.yml
conda activate rvfv-phylo
```


## Usage

For minimal pipeline options, use the ```--help``` flag e.g. 

```nextflow run main.nf --help```

To see all the options, use the ```--show_hidden_params``` flag e.g.

```nextflow run main.nf --help --show_hidden_params```

1. A typical command to perform genetic and phylogenetic analysis on then M
(medium) segment with the MP-12 strain as reference is shown as:

```
BASEDIR="${HOME}/projects/vaccine-and-circulating-strains-analysis/segments/M/complete/global"

```
```
nextflow run main.nf \
    --subworkflow rvfvcirculatingstrains \
    --fasta ${BASEDIR}/merged-sequences/RVFV-M.fasta \
    --metadata ${BASEDIR}/merged-sequences/RVFV-M.csv \
    --segment M \
    --prefix riftM \
    --start 20 \
    --end 3611 \
    --remove-duplicates false \
    --lineages ${BASEDIR}/assignment/output-dir/report/lineages.csv \
    --vaccine_reference DQ380208 \
    --outliers ./assets/RVFV-M-outliers-circulating.txt \
    --recombinants ./assets/rvfv-M-potential-recombinants.txt \
    --outdir ./output-dir/riftM \
    -work-dir ./work/riftM \
    -resume
```

2. This should be folowed by mutational profiling where the ZH548 is used a
   reference sequence.

A typical command on the analysis of the (medium) segment sequnces againts ZH548
strain is shown as:

```
nextflow run main.nf \
    --subworkflow rvfvmutationalprofiling \
    --fasta ${BASEDIR}/merged-sequences/rvfv-M.fasta \
    --metadata ${BASEDIR}/merged-sequences/rvfv-M.csv \
    --segment M \
    --prefix riftM \
    --start 20 \
    --end 3611 \
    --remove-duplicates false \
    --lineages ${BASEDIR}/assignment/output-dir/report/lineages.csv \
    --vaccine_reference NC_014396 \
    --outliers ./assets/RVFV-M-outliers-circulating.txt \
    --recombinants ./assets/rvfv-M-potential-recombinants.txt \
    --outdir ./output-dir/mutational-profiling/riftM \
    -work-dir ./work/mutational-profiling/riftM \
    -resume
```

## Method details

The pipeline offers several parameters including as highlighted:

```
Input/output options
  --subworkflow                [string]  Subworkflow type. options are 'rvfvcirculatingstrains', 'rvfvphylocontinuous' and 'rvfvmutationalprofiling'
  --fasta                      [string]  Input Fasta file containing the sequences
  --metadata                   [string]  Input comma-separated values (csv) metadata file containing the columns 'sample_name' and  'Ct.
  --recombinants               [string]  Input text file containing recombinant sequence accessions and the the column 'recombinants'
  --segment                    [string]  genomic segment of the virus. options are 'S', 'M' and 'L'
  --outdir                     [string]  The output directory where the results will be saved. [default: ./results]
  --email                      [string]  Email address for completion summary.

Alignment masking options
  --start                      [integer] start position to trim alignment (0-based index)
  --end                        [integer] end position to trim alignment (0-based index) [default: 20]

Alignment filtering options
  --outliers                   [string]  Input text file containing sequence identifiers for outlier sequences The identifiers should be in the format of the 
                                         reformatted headers as generated by REFORMAT_HEADERS process 

strain characterization options
  --lineages                   [string]  Input text file in CSV format (generatated by rvfvtyper pipeline) containing lineage information.
  --vaccine_reference          [string]  Reference accession to the vaccine strain (these accession should be in your sequence dataset)
  --grouping_column            [string]  Column to use to for grouping of sequences [default: strain_type]
  --min_freq                   [number]  Minimum percentage of sequences required to support a SNP call [default: 0.2]
  --max_freq                   [number]  Maximum percentage of sequences required to support a SNP call [default: 0.8]
  --seq_type                   [string]  Sequence type either dna or protein [default: protein]
  --snp_type                   [string]  SNP type, singleton (only a single snp per position), multiple (more than one snp per position) and conserved (snps tha 
                                         occur commonly across all the sequences) [default: singleton] 
  --group_per_lineage          [string]  Specify if you want group the stats output per lineage [default: false]

dataset filtering options
  --filter_columns             [string]  Column names (separated by space) to be used as filter to exclude sequence records with no information on the specified 
                                         columns: 'country', 'location', 'host', 'date'  [default: country date] 

Process skipping options
  --skip_modeltesting          [boolean] Skip model tesing step using modeltest-ng.

```

## Output

For ***rvfvcirculatingstrains** the typical outputs are as displayed in the tree
structure below

```
output-dir/riftM/
├── alignment
│   ├── riftM.log
│   ├── riftM_align.fasta
│   ├── riftM_dedup.fasta
│   ├── riftM_dedup.txt
│   ├── riftM_duplicated_taxa.txt
│   ├── riftM_filtered.fasta
│   ├── riftM_filtered.txt
│   └── riftM_masked.fasta
├── hyphy
│   ├── riftM.FEL.json
│   ├── riftM.FEL.log
│   ├── riftM.FUBAR.cache
│   ├── riftM.FUBAR.json
│   ├── riftM.FUBAR.log
│   ├── riftM.SLAC.json
│   └── riftM.SLAC.log
├── iqtree
│   ├── riftM.bionj
│   ├── riftM.iqtree
│   └── riftM.treefile
├── models
│   ├── riftM.model.ckp
│   ├── riftM.model.log
│   ├── riftM.model.out
│   ├── riftM.model.topos
│   └── riftM.model.tree
├── phylogeo
│   ├── riftM_dates.csv
│   ├── riftM_geocoded.pdf
│   ├── riftM_geocoded.txt
│   └── riftM_geolocations.csv
├── pipeline_info
│   ├── execution_report.html
│   ├── execution_timeline.html
│   ├── execution_trace.txt
│   └── pipeline_dag.svg
├── reformatted-sequences
│   ├── riftM.fasta
│   └── riftM.txt
├── sequences
│   ├── riftM.csv
│   └── riftM.fasta
├── strain-types
│   ├── riftM.DQ380208.all.csv
│   ├── riftM.DQ380208.per.lineage.csv
│   ├── riftM.DQ380208.snps.csv
│   ├── riftM.DQ380208.snps.pdf
│   ├── riftM.DQ380208.txt
│   ├── riftM.fasta
│   ├── riftM.strain_type.DQ380208.amino.acid.fasta
│   ├── riftM.strain_type.DQ380208.bcftools.stats.tstv.txt
│   ├── riftM.strain_type.DQ380208.mutations.per.strain.singleton.csv
│   ├── riftM.strain_type.DQ380208.parsed.bcftools.stats.csv
│   ├── riftM.strain_type.DQ380208.singleton.txt
│   ├── riftM.strain_type.DQ380208.snpSift.tstv.txt
│   ├── riftM.strain_type.DQ380208.vcf
│   ├── riftM.traits.txt
│   ├── strain_type.DQ380208.sorted.alignment.all.labels.csv
│   ├── strain_type.DQ380208.sorted.alignment.fasta
└── treetime
    ├── divergence_tree.nexus
    ├── riftM.log
    ├── riftM_ancestral_sequences.fasta
    ├── riftM_timetree.nexus
    └── trace_run.log

```
For ***rvfvmutationalprofiling** the typical outputs are as displayed in the tree
structure below

```
├── alignment
│   ├── riftM.log
│   ├── riftM.reverse.complement.fasta
│   ├── riftM_align.fasta
│   ├── riftM_dedup.fasta
│   ├── riftM_dedup.txt
│   ├── riftM_duplicated_taxa.txt
│   ├── riftM_filtered.fasta
│   ├── riftM_filtered.txt
│   └── riftM_masked.fasta
├── pipeline_info
│   ├── execution_report.html
│   ├── execution_timeline.html
│   ├── execution_trace.txt
│   └── pipeline_dag.svg
├── reformatted-sequences
│   ├── riftM.fasta
│   └── riftM.txt
├── sequence-contexts
│   ├── A3A_A3B.CT.riftM.context.txt
│   ├── A3A_A3B.GA.riftM.context.txt
│   ├── A3A_A3B.riftM.edited.sites.txt
│   ├── A3C_A3F.CT.riftM.context.txt
│   ├── A3C_A3F.GA.riftM.context.txt
│   ├── A3C_A3F.riftM.edited.sites.txt
│   ├── A3G.CT.riftM.context.txt
│   ├── A3G.GA.riftM.context.txt
│   ├── A3G.riftM.edited.sites.txt
│   ├── antelope.CT.riftM.sequence-context.txt
│   ├── antelope.GA.riftM.sequence-context.txt
│   ├── bat.CT.riftM.sequence-context.txt
│   ├── bat.GA.riftM.sequence-context.txt
│   ├── buffalo.CT.riftM.sequence-context.txt
│   ├── buffalo.GA.riftM.sequence-context.txt
│   ├── cow.CT.riftM.sequence-context.txt
│   ├── cow.GA.riftM.sequence-context.txt
│   ├── human.CT.riftM.sequence-context.txt
│   ├── human.GA.riftM.sequence-context.txt
│   ├── mosquito.CT.riftM.sequence-context.txt
│   ├── mosquito.GA.riftM.sequence-context.txt
│   ├── riftM-CT-sequence-contexts-with-metadata.csv
│   ├── riftM-GA-sequence-contexts-with-metadata.csv
│   ├── riftM-sequence-contexts-with-metadata.csv
│   ├── sheep.CT.riftM.sequence-context.txt
│   └── sheep.GA.riftM.sequence-context.txt
├── sequences
│   ├── riftM.csv
│   └── riftM.fasta
├── strain-types
│   ├── riftM.NC_014396.txt
│   ├── riftM.fasta
│   ├── riftM.strain_type.NC_014396.mutations.per.strain.singleton.csv
│   ├── riftM.strain_type.NC_014396.singleton.txt
│   ├── riftM.strain_type.NC_014396.vcf
│   ├── riftM.traits.txt
│   ├── strain_type.NC_014396.sorted.alignment.all.labels.csv
│   └── strain_type.NC_014396.sorted.alignment.fasta
└── strain-types-rev
    ├── riftM.NC_014396.txt
    ├── riftM.fasta
    ├── riftM.strain_type.NC_014396.mutations.per.strain.singleton.csv
    ├── riftM.strain_type.NC_014396.singleton.txt
    ├── riftM.strain_type.NC_014396.vcf
    ├── riftM.traits.txt
    ├── strain_type.NC_014396.sorted.alignment.all.labels.csv
    └── strain_type.NC_014396.sorted.alignment.fasta
```

3. Use the `bin/getSequenceContexts.py` to extract the sequence contexts of the
   strands from the alignments. These can be visualized using the `plots.R` to
   generate the figures.

## Credits

rvfvphylo was originally written by ajodeh-juma.

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on [Slack](https://nfcore.slack.com/channels/rvfvcirculatingstrains) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citation

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi. -->
<!-- If you use  rvfvphylo for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).  
> ReadCube: [Full Access Link](https://rdcu.be/b1GjZ)
