process EXTRACT_MRNA_SEQUENCES {
    tag "extract_mrna"
    label 'process_single'

    conda "bioconda::gffread=0.12.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gffread:0.12.7--h077b44d_5' :
        'biocontainers/gffread:0.12.7--h077b44d_6' }"

    input:
    path fasta
    path gtf

    output:
    path "mrna_sequences.fa"     , emit: fasta
    path "mrna_sequence_summary.tsv", emit: summary
    path "versions.yml"          , emit: versions

    script:
    """
    awk -F'\\t' '
    \$0 !~ /^#/ && \$9 ~ /transcript_id "/ &&
    \$9 ~ /(transcript_biotype|transcript_type|gene_biotype|gene_type|biotype) "protein_coding"/ {
        transcript_id = \$9
        sub(/.*transcript_id "/, "", transcript_id)
        sub(/".*/, "", transcript_id)
        print transcript_id
    }
    ' $gtf | sort -u > mrna_ids.txt

    awk -F'\\t' '
    FNR == NR { ids[\$1] = 1; next }
    \$0 ~ /^#/ { next }
    {
        transcript_id = \$9
        sub(/.*transcript_id "/, "", transcript_id)
        sub(/".*/, "", transcript_id)
        if (transcript_id in ids) print
    }
    ' mrna_ids.txt $gtf > mrna.gtf

    selected_transcripts=\$(wc -l < mrna_ids.txt | tr -d ' ')
    if [ "\$selected_transcripts" -lt 1 ]; then
        echo "ERROR: No protein-coding transcripts found in $gtf for FEELnc mRNA reference." >&2
        exit 1
    fi

    if [ ! -s mrna.gtf ]; then
        echo "ERROR: No GTF records found for selected protein-coding transcript IDs." >&2
        exit 1
    fi

    gffread -w mrna_sequences.fa -g $fasta mrna.gtf

    fasta_sequences=\$(grep -c '^>' mrna_sequences.fa 2>/dev/null || true)
    if [ "\${fasta_sequences:-0}" -lt 1 ]; then
        echo "ERROR: mRNA FASTA extraction produced no sequences. Check FASTA/GTF contig names and transcript annotations." >&2
        exit 1
    fi

    cat <<-END_SUMMARY > mrna_sequence_summary.tsv
    sequence_set\tselected_transcripts\tfasta_sequences
    mrna\t\${selected_transcripts}\t\${fasta_sequences:-0}
    END_SUMMARY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: \$(gffread --version 2>&1)
    END_VERSIONS
    """

    stub:
    """
    touch mrna_sequences.fa

    cat <<-END_SUMMARY > mrna_sequence_summary.tsv
    sequence_set\tselected_transcripts\tfasta_sequences
    mrna\t0\t0
    END_SUMMARY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: "0.12.7"
    END_VERSIONS
    """
}
