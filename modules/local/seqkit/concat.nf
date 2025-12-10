nextflow.preview.types = true

process SEQKIT_CONCAT {
    tag "${meta1.id}"
    label 'process_medium'

    conda "bioconda::seqkit=2.1.0"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/seqkit:2.1.0--h9ee0642_0'
        : 'quay.io/biocontainers/seqkit:2.1.0--h9ee0642_0'}"

    input:
    (meta1, fasta1) : Tuple<?, Path>
    (meta2, fasta2) : Tuple<?, Path>

    output:
    fasta    = tuple(meta1, file("*.fa"))
    versions = file("versions.yml")

    when:
    task.ext.when == null || task.ext.when

    script:
    def _args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta1.id}"
    """
    seqkit \\
        concat \\
        ${fasta1} \\
        ${fasta2} \\
        > ${prefix}_concat.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$( seqkit | sed '3!d; s/Version: //' )
    END_VERSIONS
    """
}
