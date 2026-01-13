include { STRINGTIE_STRINGTIE } from '../../../modules/nf-core/stringtie/stringtie/main'
include { STRINGTIE_MERGE      } from '../../../modules/nf-core/stringtie/merge/main'
include { SALMON_INDEX } from '../../../modules/nf-core/salmon/index/main'
include { SALMON_QUANT } from '../../../modules/nf-core/salmon/quant/main'
include { SUBREAD_FEATURECOUNTS as SUBREAD_FEATURECOUNTS_GENES } from '../../../modules/nf-core/subread/featurecounts/main'
include { SUBREAD_FEATURECOUNTS as SUBREAD_FEATURECOUNTS_TRANSCRIPTS } from '../../../modules/nf-core/subread/featurecounts/main'

workflow QUANTIFICATION {

    take:
    reads                           // channel: reads: val(meta), path(fastq)
    fasta                           // channel: fasta: val(meta), path(fasta)
    transcript_fasta                // channel: transcript_fasta: path(fasta)
    bam                             // channel: output from ALIGNMENT: val(meta), path(bam)
    reference_gtf                   // channel: reference GTF: path(gtf)

    main:
    versions = Channel.empty()

    // Initialize output channels as empty
    featurecounts_genes_out         =   Channel.empty()
    featurecounts_transcripts_out   =   Channel.empty()
    salmon_out                      =   Channel.empty()

    // Perform Salmon quantification if selected
    if (params.quantification_tool == 'featurecounts') {

        //
        // Run Stringtie
        //
        STRINGTIE_STRINGTIE(
            bam,
            reference_gtf,
        )
        versions = versions.mix(STRINGTIE_STRINGTIE.out.versions)

        // Create channel for stringtie_stringtie output
        stringtie_gtf = STRINGTIE_STRINGTIE.out.transcript_gtf.collect{it[1]}

        //
        // Merge Stringtie results
        //
        STRINGTIE_MERGE(
            stringtie_gtf,
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
    }

    // Perform Salmon quantification if selected
    if (params.quantification_tool == 'salmon') {

        // Index FASTAs
        SALMON_INDEX(
            fasta.collect{it[1]},
            transcript_fasta,
        )
        salmon_index = SALMON_INDEX.out.index

        versions = versions.mix(SALMON_INDEX.out.versions)

        // Quantification
        SALMON_QUANT(
            reads,
            salmon_index.first(),
            reference_gtf,
            transcript_fasta.first(),
            Channel.value(false),
            Channel.value(false),
        )

        versions = versions.mix(SALMON_QUANT.out.versions)

        salmon_out = SALMON_QUANT.out.results

    }

    emit:
    versions

    featurecounts_genes_out
    featurecounts_transcripts_out
    salmon_out
}
