nextflow.preview.types = true

process SUBSEQ {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/biopython_python_pyyaml:af5270f8003cfae1'
        : 'community.wave.seqera.io/library/biopython_python_pyyaml:07336da3f95d6451'}"

    input:
    (meta, fasta)        : Tuple<?, Path>
    gene                 : String
    (meta2, filter_list) : Tuple<?, Path>
    yml                  : Path

    output:
    fasta    = tuple(meta, gene, file("*.fa"))
    versions = file("versions.yml")

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def filter_arg = filter_list ? "--filter-list ${filter_list}" : ""
    def yml_arg = yml ? "--yml ${yml}" : ""
    """
    subsequence.py \\
        --sequence ${fasta} \\
        --gene ${gene} \\
        ${filter_arg} \\
        ${yml_arg} \\
        --output ${prefix}_subseq.${gene}.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
        pyyaml: \$(python -c "import yaml; print(yaml.__version__)")
    END_VERSIONS
    """
}
