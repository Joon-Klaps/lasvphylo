//modified from the original module

process IQTREE {
    tag "$alignment"
    label 'process_medium'

    conda "bioconda::iqtree=2.1.4_beta"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/iqtree:2.1.4_beta--hdcc8f71_0' :
        'quay.io/biocontainers/iqtree:2.1.4_beta--hdcc8f71_0' }"

    input:
    tuple val(meta), path(alignment)
    tuple val(meta2), path(constrain_tree)
    val constant_sites

    output:
    path "*.treefile",    emit: phylogeny
    path "*.log",         emit: log
    path "*.mldist",      emit: mldist, optional:true
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def fconst_args = constant_sites ? "-fconst $constant_sites" : ''
    def fconst_tree = constrain_tree ? "-g $constrain_tree" : ''
    def memory      = task.memory.toString().replaceAll(' ', '')
    """
    iqtree \\
        $fconst_args \\
        $args \\
        -s $alignment \\
        $fconst_tree \\
        -nt AUTO \\
        -ntmax $task.cpus \\
        -mem $memory \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iqtree: \$(echo \$(iqtree -version 2>&1) | sed 's/^IQ-TREE multicore version //;s/ .*//')
    END_VERSIONS
    """
}
