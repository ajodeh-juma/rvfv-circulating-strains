// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GET_SNPS_RECORDS {
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
    tuple val(meta), val(fasta)
    tuple val(meta), val(metadata)
    val reference
    val grouping_column
    path lineages

    output:
    tuple val(meta), path("*.csv"),   emit: mutations
    tuple val(meta), path("*.txt"),   emit: txt

    script:
    def prefix    = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    """
    get_snps_records.py \\
        --alignment $fasta \\
        --metadata $metadata \\
        --reference $reference \\
        $options.args \\
        --prefix ${prefix}.${grouping_column}.${reference} \\
        --outdir .
    """
}
