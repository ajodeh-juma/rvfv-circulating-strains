// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MODEL_TEST {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? 'bioconda::mafft=7.475' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/mafft:7.475--h779adbc_1'
    } else {
        container 'quay.io/biocontainers/mafft:7.475--h779adbc_1'
    }

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path('*.log')                           , emit: log
    tuple val(meta), path('*.out' )                          , emit: out
    tuple val(meta), path('*.topos' )                        , emit: topology
    tuple val(meta), path('*.ckp' )                          , emit: checkpoint
    tuple val(meta), path('*.tree' )                         , emit: tree
    tuple val(meta), path('*.version.txt')                   , emit: version

    script:
    def software  = getSoftwareName(task.process)
    def prefix    = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    """
    modeltest-ng \\
        --input $fasta \\
        --datatype nt \\
        -p $task.cpus \\
        --models HKY,GTR,JC,TPM1,TPM2,TPM3,TIM1,TIM2,TIM3,TVM \\
        -t ml \\
        -o ${prefix}.model 2> \\
        ${prefix}.model.log
        
    echo \$(modeltest-ng --version | grep 'ModelTest-NG' | uniq | cut -d " " -f 2) > ${software}.version.txt
    """
}