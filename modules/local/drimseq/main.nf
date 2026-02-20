process DRIMSEQ {
    label 'process_medium'

    conda "bioconda::bioconductor-drimseq=1.28.0 bioconda::bioconductor-tximport=1.28.0 bioconda::bioconductor-genomicfeatures=1.52.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-drimseq:1.34.0--r44hdfd78af_0' :
        'biocontainers/bioconductor-drimseq:1.34.0--r44hdfd78af_0' }"

    input:
    path count_files
    path gtf
    path design_file

    output:
    path "*.pdf"                   , optional:true, emit: plots
    path "*.csv"                   , emit: tables
    path "*.rds"                   , emit: rds
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'drimseq.R'
}
