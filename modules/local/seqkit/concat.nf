process SEQKIT_CONCAT{
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::seqkit=2.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.1.0--h9ee0642_0':
        'quay.io/biocontainers/seqkit:2.1.0--h9ee0642_0' }"

    input:
    path fasta1
    path fasta2

    output:
    path "*.fa"         , emit: fasta
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    sed 's/>_R_/>/g' ${fasta1} > fasta1.fas
    sed 's/>_R_/>/g' ${fasta2} > fasta2.fas
    seqkit \\
        concat \\
        -f fasta1.fas \\
        fasta2.fas \\
        > ${prefix}_concat.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$( seqkit | sed '3!d; s/Version: //' )
    END_VERSIONS
    """
}
