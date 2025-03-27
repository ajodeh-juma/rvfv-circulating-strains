// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SNIPIT {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::snipit=1.1.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/snipit:1.1.2--pyhdfd78af_0"
    } else {
        container "quay.io/biocontainers/snipit:1.1.2--pyhdfd78af_0"
    }

    input:
    tuple val(meta), path(fasta)
    val reference
    tuple val(meta), path(metadata)
    tuple val(meta), path(snps)
    

    output:
    tuple val(meta), path("*.pdf")    , emit: pdf
    tuple val(meta), path("*.csv")    , emit: csv

    script:
    def software          = getSoftwareName(task.process)
    def memory            = task.memory.toString().replaceAll(' ', '')
    def prefix            = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
        
    """
    snipit $fasta \\
        -r $reference \\
        -d . \\
        -o ${prefix}.${reference}.snps \\
        -s \\
        --labels $metadata \\
        --exclude-ambig-pos \\
        --l-header 'name,label' \\
        --include-positions \$(< ${snps}) \\
        --colour-palette purine-pyrimidine \\
        --solid-background \\
        --format pdf \\
        --height 32 \\
        --width 36 \\
        --write-snps \\
        --size-option scale

    mv snps.csv ${prefix}.${reference}.snps.csv
    """
    
}
