// Modified version from the original MAFFT nf-core module

process MAFFT {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::mafft=7.508"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mafft:7.508--hec16e2b_0':
        'quay.io/biocontainers/mafft:7.508--hec16e2b_0' }"

    input:
    tuple val(meta), path(addsequences), path(alignment), path(tree)
    path fasta
    val(onlynewseq) // boolean


    output:
    tuple val(meta), path("*.fas"), path(alignment), path(tree)         , emit: addseq
    path ("*.fas")                                                      , emit: fasta
    path "versions.yml"                                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def fasta = fasta ?: "${alignment}" // if external alginment not given use the one in the map
    def prefix = task.ext.prefix ?: "${meta.id}"
    def add = addsequences ? "--add ${addsequences}" : ''
    def onlynewseq = onlynewseq ?: 'FALSE'
    """
    mafft \\
        --thread ${task.cpus} \\
        ${args} \\
        ${add} \\
        ${fasta} \\
        > ${prefix}.fas

    if ${onlynewseq}; then 
        selectSeqs.sh ${prefix}.fas ${add}
    fi 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mafft: \$(mafft --version 2>&1 | sed 's/^v//' | sed 's/ (.*)//')
    END_VERSIONS
    """
}
