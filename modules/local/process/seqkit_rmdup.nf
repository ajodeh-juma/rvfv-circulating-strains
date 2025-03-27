// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SEQKIT_RMDUP {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? 'bioconda::seqkit=2.4.0' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/seqkit:2.4.0--h527b516_0'
    } else {
        container 'quay.io/biocontainers/seqkit:2.4.0--h527b516_0'
    }

    input:
    tuple val(meta), path(alignment)

    output:
    tuple val(meta), path('*.fasta')                           , emit: fasta
    tuple val(meta), path('*.txt' )                            , emit: txt
    // tuple val(meta), path('*.version.txt')                     , emit: version

    script:
    def software  = getSoftwareName(task.process)
    def prefix    = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    """
    seqkit rmdup -s -D  ${prefix}_dedup.txt < $alignment > ${prefix}_dedup.fasta

    if [ ! -f ${prefix}_dedup.txt ]; then
        touch ${prefix}_dedup.txt
    fi
    
    """
    // echo \$(mafft --version 2>&1) | sed -e 's/^*.v\$//' > ${software}.version.txt
}