// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SNPSITES {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::snp-sites=2.5.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/snp-sites:2.5.1"
    } else {
        container "quay.io/biocontainers/snp-sites:2.5.1"
    }

    input:
    tuple val(meta), path(fasta)
    val grouping_column
    val reference

    output:
    tuple val(meta), path("*.vcf"), emit: vcf

    script:
    def prefix    = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    """
    snp-sites -vc -o ${prefix}.${grouping_column}.${reference}.vcf $fasta
    """
}
