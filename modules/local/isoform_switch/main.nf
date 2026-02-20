process ISOFORM_SWITCH {
    label 'process_medium'

    conda "bioconda::bioconductor-isoformswitchanalyzer=2.0.0 bioconda::bioconductor-drimseq=1.28.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-isoformswitchanalyzer:2.2.0--r43ha9d7317_0' :
        'biocontainers/bioconductor-isoformswitchanalyzer:2.6.0--r44h3df3fcb_0' }"

    input:
    path count_files
    path tpm_files
    path fasta
    path gtf
    path design_file
    path gene_results_file
    path transcript_results_file

    output:
    path "isoform_switch_results.rds", emit: rds
    path "*.pdf"                     , optional:true, emit: plots
    path "*.csv"                     , optional:true, emit: tables
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'isoform_switch.R'
}
