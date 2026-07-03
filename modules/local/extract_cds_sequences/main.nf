process EXTRACT_CDS_SEQUENCES {
    tag "extract_cds"
    label 'process_single'

    conda "bioconda::gffread=0.12.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gffread:0.12.7--h077b44d_5' :
        'biocontainers/gffread:0.12.7--h077b44d_6' }"

    input:
    path fasta
    path gtf

    output:
    path "cds_sequences.fa"             , emit: fasta
    path "cds_sequence_summary.tsv"     , emit: summary
    path "versions.yml"                 , emit: versions

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
    ' $gtf | sort -u > cds_ids.txt

    awk -F'\\t' '
    FNR == NR { ids[\$1] = 1; next }
    \$0 ~ /^#/ { next }
    {
        transcript_id = \$9
        sub(/.*transcript_id "/, "", transcript_id)
        sub(/".*/, "", transcript_id)
        if (transcript_id in ids) print
    }
    ' cds_ids.txt $gtf > cds.gtf

    selected_transcripts=\$(wc -l < cds_ids.txt | tr -d ' ')
    if [ "\$selected_transcripts" -lt 1 ]; then
        echo "ERROR: No protein-coding transcripts found in $gtf for CPAT coding training." >&2
        exit 1
    fi

    if [ ! -s cds.gtf ]; then
        echo "ERROR: No GTF records found for selected protein-coding transcript IDs." >&2
        exit 1
    fi

    gffread -x cds_sequences.fa -g $fasta cds.gtf

    fasta_sequences=\$(grep -c '^>' cds_sequences.fa 2>/dev/null || true)
    if [ "\${fasta_sequences:-0}" -lt 1 ]; then
        echo "ERROR: CDS FASTA extraction produced no sequences. Check FASTA/GTF contig names and CDS annotations." >&2
        exit 1
    fi

    cat <<-END_SUMMARY > cds_sequence_summary.tsv
    sequence_set\tselected_transcripts\tfasta_sequences
    cds\t\${selected_transcripts}\t\${fasta_sequences:-0}
    END_SUMMARY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: \$(gffread --version 2>&1)
    END_VERSIONS
    """

    stub:
    """
    touch cds_sequences.fa

    cat <<-END_SUMMARY > cds_sequence_summary.tsv
    sequence_set\tselected_transcripts\tfasta_sequences
    cds\t0\t0
    END_SUMMARY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: "0.12.7"
    END_VERSIONS
    """
}
