// Modified version from the original MAFFT nf-core module

process MAFFT {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::mafft=7.508"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mafft:7.508--hec16e2b_0':
        'quay.io/biocontainers/mafft:7.508--hec16e2b_0' }"

    input:
    tuple val(meta), path(addsequences), path(alignment)
    val gene


    output:
    tuple val(meta), path("*.fas")  , emit: fasta
    path ("*.txt")                  , emit: pattern
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def fasta = alignment 
    def prefix = task.ext.prefix ?: "${meta.id}"
    def add = addsequences ? "--add ${addsequences}" : ''
    def gene = gene? ".${gene}" : ''
    
    """
    mafft \\
        --thread ${task.cpus} \\
        ${args} \\
        ${add} \\
        ${fasta}  \\
        > tmp

    sed "s/>_R_/>/g" tmp > ${prefix}${gene}.fas 
    grep ">" ${addsequences} | sed "s/>//g" > ${prefix}.patterns.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mafft: \$(mafft --version 2>&1 | sed 's/^v//' | sed 's/ (.*)//')
    END_VERSIONS
    """
}
