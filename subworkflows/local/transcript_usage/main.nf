//
// Subworkflow for Transcript Usage Analysis
//

include { DRIMSEQ        } from '../../../modules/local/drimseq/main'
include { DEXSEQ         } from '../../../modules/local/dexseq/main'
include { ISOFORM_SWITCH } from '../../../modules/local/isoform_switch/main'

workflow TRANSCRIPT_USAGE {

    take:
    transcript_counts // channel: [ meta, path(counts_transcript) ]
    transcript_tpm    // channel: [ meta, path(tpm_transcript) ]
    fasta             // channel: path(fasta)
    gtf               // channel: path(gtf)
    samples           // channel: [ meta, reads ]

    main:

    ch_versions = Channel.empty()

    //
    // Create design file
    //
    samples
        .map { meta, reads ->
            "${meta.id},${meta.condition}"
        }
        .collectFile(name: 'design.csv', newLine: true, sort: true, seed: "sampleID,condition")
        .set { ch_design_file }

    //
    // Collect all quantification files
    //
    transcript_counts
        .map{ it[1] }
        .collect()
        .set { ch_count_files }

    transcript_tpm
        .map{ it[1] }
        .collect()
        .set { ch_tpm_files }

    //
    // Run Statistical Testing (DRIMSeq or DEXSeq)
    //
    ch_gene_results = Channel.empty()
    ch_transcript_results = Channel.empty()
    ch_drimseq_tables = Channel.empty()
    ch_drimseq_plots = Channel.empty()

    if (params.dtu_tool == 'dexseq') {
        DEXSEQ(
            ch_count_files,
            gtf,
            ch_design_file
        )
        ch_gene_results = DEXSEQ.out.tables.map{ it -> [it.find{f -> f.name == "dexseq_gene_results.csv"}] }
        ch_transcript_results = DEXSEQ.out.tables.map{ it -> [it.find{f -> f.name == "dexseq_transcript_results.csv"}] }
        ch_versions = ch_versions.mix(DEXSEQ.out.versions)
    } else {
        // Default to DRIMSeq
        DRIMSEQ(
            ch_count_files,
            gtf,
            ch_design_file
        )
        ch_gene_results = DRIMSEQ.out.tables.map{ it -> [it.find{f -> f.name == "drimseq_gene_results.csv"}] }
        ch_transcript_results = DRIMSEQ.out.tables.map{ it -> [it.find{f -> f.name == "drimseq_transcript_results.csv"}] }
        ch_drimseq_tables = DRIMSEQ.out.tables
        ch_drimseq_plots = DRIMSEQ.out.plots
        ch_versions = ch_versions.mix(DRIMSEQ.out.versions)
    }

    //
    // Run IsoformSwitchAnalyzeR
    //
    ISOFORM_SWITCH(
        ch_count_files,
        ch_tpm_files,
        fasta,
        gtf,
        ch_design_file,
        ch_gene_results,
        ch_transcript_results
    )
    ch_versions = ch_versions.mix(ISOFORM_SWITCH.out.versions)

    emit:
    drimseq_tables      = ch_drimseq_tables
    drimseq_plots       = ch_drimseq_plots
    isoform_switch_plots = ISOFORM_SWITCH.out.plots
    isoform_switch_tables = ISOFORM_SWITCH.out.tables
    versions            = ch_versions
}
