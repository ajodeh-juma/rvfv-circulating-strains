// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process HYPHY_SLAC {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? 'bioconda::hyphy=2.5.38' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/hyphy:2.3.58--ha272830_0'
    } else {
        container 'quay.io/biocontainers/hyphy:2.3.58--ha272830_0'
    }

    input:
    tuple val(meta), path(fasta)
    tuple val(meta), path(treefile)


    output:
    tuple val(meta), path('*.json')                            , emit: json
    tuple val(meta), path('*.log')                             , emit: log
    tuple val(meta), path('*.version.txt')                     , emit: version

    script:
    def software  = getSoftwareName(task.process)
    def prefix    = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    """
    hyphy slac \\
      CPU=$task.cpus \\
      -type Nucleotide \\
      --alignment $fasta \\
      --tree $treefile \\
      --rv Gamma \\
      --rate-classes 4 \\
      --mode Faster \\
      --output ${prefix}.SLAC.json \\
      2>&1 > ${prefix}.SLAC.log
      

    echo \$(hyphy --version | cut -d " " -f2) > ${software}.version.txt
    """
}