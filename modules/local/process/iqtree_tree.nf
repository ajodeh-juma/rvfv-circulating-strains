// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process IQTREE_TREE {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::iqtree=2.0.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/iqtree:2.0.3--h176a8bc_1"
    } else {
        container "quay.io/biocontainers/iqtree:2.0.3--h176a8bc_1"
    }

    input:
    tuple val(meta), path(alignment)
    tuple val(meta), path(models)

    output:
    tuple val(meta), path("*.treefile")    , emit: treefile
    tuple val(meta), path("*.bionj")       , emit: consensus
    tuple val(meta), path("*.iqtree")      , emit: iqtree
    path "*.version.txt"                   , emit: version

    script:
    def software          = getSoftwareName(task.process)
    def memory            = task.memory.toString().replaceAll(' ', '')
    def prefix            = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
        
    """
    MODEL=`grep 'iqtree' ${models} | cut -d " " -f8 | uniq -c | sort -nr | head -n 1 | awk -F" " '{ print \$2 }'`

    iqtree \\
        -s $alignment \\
        -T $task.cpus \\
        -m \$MODEL \\
        -mem $memory \\
        $options.args \\
        --prefix ${prefix}

    echo \$(iqtree -version 2>&1) | sed 's/^IQ-TREE multicore version \\([0-9\\.]*\\) .*\$/\\1/' > ${software}.version.txt
    """
    
}
