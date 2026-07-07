include { STRINGTIE2                                                    } from '../../../modules/local/stringtie2/main'
include { STRINGTIE_MERGE                                               } from '../../../modules/nf-core/stringtie/merge/main'
include { SUBREAD_FEATURECOUNTS as SUBREAD_FEATURECOUNTS_GENES          } from '../../../modules/nf-core/subread/featurecounts/main'
include { SUBREAD_FEATURECOUNTS as SUBREAD_FEATURECOUNTS_TRANSCRIPTS    } from '../../../modules/nf-core/subread/featurecounts/main'

workflow STRINGTIE_FEATURECOUNTS {

    take:
    fasta                           // channel: reference FASTA: path(fasta)
    bam                             // channel: output from ALIGNMENT: val(meta), path(bam)
    reference_gtf                   // channel: reference GTF: path(gtf)

    main:
    versions = Channel.empty()

    // Initialize output channels as empty
    featurecounts_genes_out         =   Channel.empty()
    featurecounts_transcripts_out   =   Channel.empty()

    //
    // Run Stringtie
    //
    fasta_file = fasta.map { it[1] }

    stringtie2_input = bam
        .combine(fasta_file)
        .combine(reference_gtf)
        .map { meta, bam_file,fasta_f, gtf ->
            tuple(meta, fasta_f, gtf, bam_file)
        }

    STRINGTIE2(
        stringtie2_input
    )
    versions = versions.mix(STRINGTIE2.out.versions)

    // Create channel for StringTie per-sample output
    stringtie_gtf = STRINGTIE2.out.stringtie_gtf
    stringtie_gtf_for_merge = stringtie_gtf.collect{it}

    //
    // Merge Stringtie results
    //
    STRINGTIE_MERGE(
        stringtie_gtf_for_merge,
        reference_gtf,
    )
    versions = versions.mix(STRINGTIE_MERGE.out.versions)

    // Create channel for stringtie_merge output
    merged_gtf = STRINGTIE_MERGE.out.gtf

    // Create channel for featurecounts input (gene level)
    Channel.value([id:'gene', single_end:true])
        .combine(bam.collect{ it[1]}.map { bams -> [bams] } )
        .combine(merged_gtf)
        .set{ featurecounts_input_genes }

    // Run quantification on gene level
    SUBREAD_FEATURECOUNTS_GENES(
        featurecounts_input_genes,
    )
    versions = versions.mix(SUBREAD_FEATURECOUNTS_GENES.out.versions)

    // Create channel for featurecounts input (transcript level)
    Channel.value([id:'transcript', single_end:true])
        .combine(bam.collect{ it[1]}.map { bams -> [bams] } )
        .combine(merged_gtf)
        .set{ featurecounts_input_transcripts }

    // Run quantification on transcript level
    SUBREAD_FEATURECOUNTS_TRANSCRIPTS(
        featurecounts_input_transcripts,
    )

    versions = versions.mix(SUBREAD_FEATURECOUNTS_TRANSCRIPTS.out.versions)
    featurecounts_genes_out = SUBREAD_FEATURECOUNTS_GENES.out.counts
    featurecounts_transcripts_out = SUBREAD_FEATURECOUNTS_TRANSCRIPTS.out.counts

    emit:
    versions

    stringtie_gtf
    merged_gtf
    featurecounts_genes_out
    featurecounts_transcripts_out
}
