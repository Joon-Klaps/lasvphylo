nextflow.preview.types = true

nextflow.preview.types = true

process MAFFT_ALIGN {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/mulled-v2-12eba4a074f913c639117640936668f5a6a01da6:425707898cf4f85051b77848be253b88f1d2298a-0'
        : 'biocontainers/mulled-v2-12eba4a074f913c639117640936668f5a6a01da6:425707898cf4f85051b77848be253b88f1d2298a-0'}"

    input:
    (meta, fasta)         : Tuple<Map, Path>
    (_meta2, add)          : Tuple<Map, Path>
    (_meta3, addfragments) : Tuple<Map, Path>
    (_meta4, addfull)      : Tuple<Map, Path>
    (_meta5, addprofile)   : Tuple<Map, Path>
    (_meta6, addlong)      : Tuple<Map, Path>
    compress              : Boolean

    output:
    fas      = tuple(meta, file("*.fas{.gz,}"))
    pattern  = tuple(meta, file("*.patterns.txt"))
    versions = file("versions.yml")

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def command_add = add ? "--add <(unpigz -cdf ${add})" : ''
    def command_addfragments = addfragments ? "--addfragments <(unpigz -cdf ${addfragments})" : ''
    def command_addfull = addfull ? "--addfull <(unpigz -cdf ${addfull})" : ''
    def command_addprofile = addprofile ? "--addprofile <(unpigz -cdf ${addprofile})" : ''
    def command_addlong = addlong ? "--addlong <(unpigz -cdf ${addlong})" : ''
    def _write_output = compress ? " | pigz -cp ${task.cpus} > ${prefix}.fas.gz" : "> ${prefix}.fas"
    // this will not preserve MAFFTs return value, but mafft crashes when it receives a process substitution
    if ("${fasta}" == "${prefix}.fas") {
        error("Input and output names are the same, set prefix in module configuration to disambiguate!")
    }
    """
    mafft \\
        --thread ${task.cpus} \\
        ${command_add} \\
        ${command_addfragments} \\
        ${command_addfull} \\
        ${command_addprofile} \\
        ${command_addlong} \\
        ${args} \\
        ${fasta} \\
        > tmp

    sed "s/>_R_/>/g" tmp > ${prefix}.fas
    grep ">" ${add} | sed "s/>//g" > ${prefix}.patterns.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mafft: \$(mafft --version 2>&1 | sed 's/^v//' | sed 's/ (.*)//')
        pigz: \$(echo \$(pigz --version 2>&1) | sed 's/^.*pigz\\w*//' ))
    END_VERSIONS
    """

    stub:
    def _args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def _command_add = add ? "--add ${add}" : ''
    def _command_addfragments = addfragments ? "--addfragments ${addfragments}" : ''
    def _command_addfull = addfull ? "--addfull ${addfull}" : ''
    def _command_addprofile = addprofile ? "--addprofile ${addprofile}" : ''
    def _command_addlong = addlong ? "--addlong ${addlong}" : ''
    if ("${fasta}" == "${prefix}.fas") {
        error("Input and output names are the same, set prefix in module configuration to disambiguate!")
    }
    """
    touch ${prefix}.fas${compress ? '.gz' : ''}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mafft: \$(mafft --version 2>&1 | sed 's/^v//' | sed 's/ (.*)//')
        pigz: \$(echo \$(pigz --version 2>&1) | sed 's/^.*pigz\\w*//' ))
    END_VERSIONS
    """
}
