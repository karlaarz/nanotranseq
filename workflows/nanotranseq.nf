/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { RAW_READS_QC                    } from '../subworkflows/local/raw_read_qc/main'
include { MULTIQC                         } from '../modules/nf-core/multiqc/main'
include { DIRECT_RNA_QC                   } from '../subworkflows/local/direct_rna_qc/main'
include { ALIGNMENT                       } from '../subworkflows/local/alignment/main'
include { BEDTOOLS_BIGWIG                 } from '../subworkflows/local/bedtools_bigwig/main'
include { STRINGTIE_FEATURECOUNTS         } from '../subworkflows/local/stringtie_featurecounts/main'
include { PSEUDOALIGNMENT                 } from '../subworkflows/local/pseudoalignment/main'
include { TRANSCRIPT_USAGE                } from '../subworkflows/local/transcript_usage/main'
include { paramsSummaryMap                } from 'plugin/nf-schema'
include { paramsSummaryMultiqc            } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML          } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText          } from '../subworkflows/local/utils_nfcore_nanotranseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NANOTRANSEQ {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    ch_fasta       // channel: fasta read in from --fasta
    ch_gtf          // channel: gtf file read in from --gtf
    ch_gene_id      // channel: attributed gene ID in the GTF file
    ch_gene_attributes     // channel: extra gene attributes in the GTF file
    ch_transcript_fasta     // channel: transcript fasta file read in from --transcript_fasta
    ch_direct_rna // channel: direct_rna read in from --direct_rna
    ch_minimap2_index   // channel: index read in from --minimap2_index

    main:

    ch_versions = channel.empty()

    //
    // MODULE: Run QC on raw reads
    //
    RAW_READS_QC (
        ch_samplesheet
    )
    ch_versions = ch_versions.mix(RAW_READS_QC.out.versions)

    //
    // Run Chopper if cDNA sequencing was performed
    if (!params.direct_rna) {

        DIRECT_RNA_QC (
            ch_samplesheet,
        )

        ch_versions = ch_versions.mix(DIRECT_RNA_QC.out.versions)

    }

    // If cDNA was performed, use CHOPPER's output as reads. If not, use raw data
    ch_reads = params.direct_rna ? ch_samplesheet : DIRECT_RNA_QC.out.reads

    //
    // Run alignment if either `featurecounts` or `both` is selected as quantification tool
    //
    if (params.quantification_tool == 'featurecounts' || params.quantification_tool == 'both') {

        // Run alignment with Minimap2
        ALIGNMENT(
            ch_reads,
            ch_fasta,
            ch_minimap2_index
        )
        ch_versions = ch_versions.mix(ALIGNMENT.out.versions)

        // Generate BigWig files for visualisation in genome browsers
        BEDTOOLS_BIGWIG(
            ch_fasta,
            ALIGNMENT.out.minimap2_bam
        )

        // Assemble and quantify
        STRINGTIE_FEATURECOUNTS(
            ch_fasta,
            ALIGNMENT.out.minimap2_bam,
            ch_gtf,
        )
        ch_versions = ch_versions.mix(STRINGTIE_FEATURECOUNTS.out.versions)

    }

    //
    // Run pseudoalignment if either `salmon` or `both` is selected as quantification tool
    //
    if (params.quantification_tool == 'salmon' || params.quantification_tool == 'both') {

        PSEUDOALIGNMENT(
            ch_reads,
            ch_fasta,
            ch_transcript_fasta,
            ch_gtf,
            ch_gene_id,
            ch_gene_attributes,
        )
        ch_versions = ch_versions.mix(PSEUDOALIGNMENT.out.versions)

        //
        // SUBWORKFLOW: Transcript Usage
        //
        TRANSCRIPT_USAGE(
            PSEUDOALIGNMENT.out.counts_transcript,
            PSEUDOALIGNMENT.out.tpm_transcript,
            ch_fasta.map{ it[1] },
            ch_gtf,
            ch_reads
        )
        ch_versions = ch_versions.mix(TRANSCRIPT_USAGE.out.versions)

    }

    //
    // Collate and save software versions
    //
    def topic_versions = Channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'nanotranseq_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    emit:
    multiqc_report = RAW_READS_QC.out.multiqc_report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
