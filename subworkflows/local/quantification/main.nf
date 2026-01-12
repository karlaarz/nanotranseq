include { STRINGTIE_STRINGTIE } from '../../../modules/nf-core/stringtie/stringtie/main'
include { STRINGTIE_MERGE      } from '../../../modules/nf-core/stringtie/merge/main'
include { SALMON_QUANT } from '../../../modules/nf-core/salmon/quant/main'
include { SUBREAD_FEATURECOUNTS as SUBREAD_FEATURECOUNTS_GENES } from '../../../modules/nf-core/subread/featurecounts/main'
include { SUBREAD_FEATURECOUNTS as SUBREAD_FEATURECOUNTS_TRANSCRIPTS } from '../../../modules/nf-core/subread/featurecounts/main'

workflow QUANTIFICATION {

    take:
    bam                             // channel: output from ALIGNMENT: val(meta), path(bam)
    reference_gtf                   // channel: reference GTF: path(gtf)

    main:
    versions = Channel.empty()
    
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

    if (params.quantification_tool == 'featurecounts') {

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

    }
    
    emit:
    versions
}