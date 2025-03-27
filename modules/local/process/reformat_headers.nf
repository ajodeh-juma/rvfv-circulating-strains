// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process REFORMAT_HEADERS {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "conda-forge::python=3.6.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/python:3.6.1"
    } else {
        container "quay.io/biocontainers/python:3.6.1"
    }

    input:
    tuple val(meta), val(files)

    output:
    tuple val(meta), path("*.fasta"), emit: fasta
    tuple val(meta), path("*.txt"),   emit: txt

    script:
    def prefix    = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    """
    reformat_headers.py \\
        --fasta ${files[0]} \\
        --metadata ${files[1]} \\
        --coordinates false \\
        --sep comma \\
        --prefix ${prefix} \\
        --headers host country date \\
        --outDir .
    """
}
