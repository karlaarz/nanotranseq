process FEELNC_CODPOT_PARSE {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.79' :
        'biocontainers/biopython:1.79' }"

    input:
    tuple val(meta), path(feelnc_rf_txt)      // *_RF.txt codpot_full
    tuple val(meta2), path(candidate_fasta)    // exon_filtered.fa

    output:
    tuple val(meta), path("*.feelnc.tsv") , emit: feelnc_results
    tuple val(meta), path("*.lncRNA.fa")  , emit: lncrna_fasta
    tuple val(meta), path("*.mRNA.fa")    , emit: mrna_fasta
    path "versions.yml"                   , emit: versions

    script:
    template 'parse_feelnc_output.py'

    stub:
    prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}.feelnc.tsv
    touch ${prefix}.lncRNA.fa
    touch ${prefix}.mRNA.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "stub"
    END_VERSIONS
    """
}
