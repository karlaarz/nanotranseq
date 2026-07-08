process EXTRACT_LNCRNA_SEQUENCES {
    tag "extract_lncrna"
    label 'process_single'

    conda "bioconda::gffread=0.12.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gffread:0.12.7--h077b44d_5' :
        'biocontainers/gffread:0.12.7--h077b44d_6' }"

    input:
    path fasta
    path gtf
    val lncrna_biotypes

    output:
    path "lncrna_sequences.fa"              , emit: fasta
    path "lncrna_sequence_summary.tsv"      , emit: summary
    path "versions.yml"                     , emit: versions

    script:
    """
    biotype_pattern=\$(printf '%s\\n' "$lncrna_biotypes" | tr ',' '\\n' | sed 's/^ *//;s/ *\$//' | awk 'NF' | paste -sd '|' -)
    if [ -z "\$biotype_pattern" ]; then
        echo "ERROR: --lncrna_biotypes is empty. Provide at least one noncoding biotype for CPAT negative training." >&2
        exit 1
    fi

    awk -v pattern="\$biotype_pattern" -F'\\t' '
    \$0 !~ /^#/ && \$9 ~ /transcript_id "/ &&
    \$9 ~ "(transcript_biotype|transcript_type|gene_biotype|gene_type|biotype) \\"(" pattern ")\\"" {
        transcript_id = \$9
        sub(/.*transcript_id "/, "", transcript_id)
        sub(/".*/, "", transcript_id)
        print transcript_id
    }
    ' $gtf | sort -u > lncrna_ids.txt

    awk -F'\\t' '
    FNR == NR { ids[\$1] = 1; next }
    \$0 ~ /^#/ { next }
    {
        transcript_id = \$9
        sub(/.*transcript_id "/, "", transcript_id)
        sub(/".*/, "", transcript_id)
        if (transcript_id in ids && \$3 != "CDS" && \$3 != "start_codon" && \$3 != "stop_codon") print
    }
    ' lncrna_ids.txt $gtf > lncrna.gtf

    selected_transcripts=\$(wc -l < lncrna_ids.txt | tr -d ' ')
    if [ "\$selected_transcripts" -lt 1 ]; then
        echo "ERROR: No lncRNA/noncoding transcripts found in $gtf using biotypes: $lncrna_biotypes" >&2
        exit 1
    fi

    if [ ! -s lncrna.gtf ]; then
        echo "ERROR: No GTF records found for the selected lncRNA/noncoding transcript IDs." >&2
        exit 1
    fi

    gffread -w lncrna_sequences.fa -g $fasta lncrna.gtf

    fasta_sequences=\$(grep -c '^>' lncrna_sequences.fa 2>/dev/null || true)
    if [ "\${fasta_sequences:-0}" -lt 1 ]; then
        echo "ERROR: lncRNA/noncoding FASTA extraction produced no sequences. Check FASTA/GTF contig names and biotypes." >&2
        exit 1
    fi

    cat <<-END_SUMMARY > lncrna_sequence_summary.tsv
    sequence_set\tselected_transcripts\tfasta_sequences
    lncrna\t\${selected_transcripts}\t\${fasta_sequences:-0}
    END_SUMMARY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: \$(gffread --version 2>&1)
    END_VERSIONS
    """

    stub:
    """
    touch lncrna_sequences.fa

    cat <<-END_SUMMARY > lncrna_sequence_summary.tsv
    sequence_set\tselected_transcripts\tfasta_sequences
    lncrna\t0\t0
    END_SUMMARY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: "0.12.7"
    END_VERSIONS
    """
}
