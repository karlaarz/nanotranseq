process DEXSEQ {
    label 'process_medium'

    conda "bioconda::bioconductor-dexseq=1.52.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-dexseq:1.48.0--r43hdfd78af_2' :
        'biocontainers/bioconductor-dexseq:1.52.0--r44hdfd78af_0' }"

    input:
    path count_files
    path gtf
    path design_file

    output:
    path "*.csv"                   , optional:true, emit: tables
    path "*.pdf"                   , optional:true, emit: plots
    path "*.rds"                   , emit: rds
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'dexseq.R'
}
