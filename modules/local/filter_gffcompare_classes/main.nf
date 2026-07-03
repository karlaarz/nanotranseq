process FILTER_GFFCOMPARE_CLASSES {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.79' :
        'biocontainers/biopython:1.79' }"

    input:
    tuple val(meta), path(gtf)
    tuple val(meta2), path(tmap)
    val class_codes

    output:
    tuple val(meta), path("*.novel.gtf")    , emit: filtered_gtf
    tuple val(meta), path("*.novel.tmap")   , emit: filtered_tmap
    path "*.class_summary.tsv"              , emit: class_summary
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'filter_gffcompare_classes.py'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.novel.gtf
    touch ${prefix}.novel.tmap
    touch ${prefix}.class_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "stub"
    END_VERSIONS
    """
}
