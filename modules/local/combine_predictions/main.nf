process COMBINE_PREDICTIONS {
    tag "combine"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.79' :
        'biocontainers/biopython:1.79' }"

    input:
    tuple val(meta1), path(cpat_results)
    tuple val(meta2), path(feelnc_results)
    tuple val(meta3), path(plek_results)
    tuple val(meta4), path(gtf)
    tuple val(meta5), path(fasta)
    path tmap

    output:
    tuple val(meta4), path("*.final_protein_coding.gtf") , emit: protein_coding_gtf
    tuple val(meta5), path("*.final_protein_coding.fa")  , emit: protein_coding_fasta
    path "*.prediction_summary.tsv"                      , emit: summary
    path "*.prediction_summary_protein_coding.tsv"       , emit: protein_coding_pred_summary
    path "*.coding_potential_report.txt"                 , emit: report
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'combine_lncrna_predictions.py'

    stub:
    def prefix = task.ext.prefix ?: "${meta4.id}"
    """
    touch ${prefix}.final_protein_coding.gtf
    touch ${prefix}.final_protein_coding.fa
    touch ${prefix}.prediction_summary.tsv
    touch ${prefix}.prediction_summary_protein_coding.tsv
    touch ${prefix}.coding_potential_report.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "stub"
    END_VERSIONS
    """
}
