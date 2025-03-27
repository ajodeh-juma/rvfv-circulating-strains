// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SUBSET_BY_COUNTRY {
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
    tuple val(meta), val(prefix)
    val fasta
    val metadata

    output:
    tuple val(meta), path("*.fasta"), emit: fasta
    tuple val(meta), path("*.csv"),   emit: csv

    script:
    //def prefix    = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    """
    subset_by_country.py \\
        --fasta $fasta \\
        --csv $metadata \\
        --countries "Burundi" "Burkina Faso" "Central African Republic" "Egypt" "Guinea" "Kenya" "Madagascar" "Mauritania" "Mayotte" "Namibia" "Rwanda" "Sudan" "Senegal" "Tanzania" "Uganda" "South Africa" "Zimbabwe" \\
        --prefix ${prefix} \\
        --outDir .
    """
}
