// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process TREETIME_TREE {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::treetime=0.9.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/treetime:0.9.3--pyh5e36f6f_0"
    } else {
        container "quay.io/biocontainers/treetime:0.9.3--pyh5e36f6f_0"
    }

    input:
    tuple val(meta), path(tree)
    tuple val(meta), path(alignment)
    tuple val(meta), path(dates)

    output:
    tuple val(meta), path("*.nexus")       , emit: nexus
    tuple val(meta), path('*.fasta')       , emit: ancestral
    tuple val(meta), path('*.log')         , emit: log
    path "*.version.txt"                   , emit: version

    script:
    def software          = getSoftwareName(task.process)
    def prefix            = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
        
    """
    treetime \\
        --tree $tree \\
        --aln $alignment \\
        --dates $dates \\
        --outdir . \\
        2>&1 > ${prefix}.log

    mv timetree.nexus ${prefix}_timetree.nexus
    mv ancestral_sequences.fasta ${prefix}_ancestral_sequences.fasta
    echo \$(treetime -version 2>&1 | cut -d " " -f2) > ${software}.version.txt
    """
    
}
