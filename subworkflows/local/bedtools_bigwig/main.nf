include { SAMTOOLS_FAIDX                } from '../../../modules/nf-core/samtools/faidx/main' 
include { BEDTOOLS_GENOMECOV            } from '../../../modules/nf-core/bedtools/genomecov/main'
include { UCSC_BEDGRAPHTOBIGWIG         } from '../../../modules/nf-core/ucsc/bedgraphtobigwig/main'

workflow BEDTOOLS_BIGWIG {
    take:
    ch_fasta
    ch_sorted_bam

    main:

    versions = channel.empty()

    /*
     * Get chromosome sizes
     */
    SAMTOOLS_FAIDX(
        ch_fasta.map { meta, fasta -> tuple(meta, fasta, []) },
        true
    )
    ch_chr_sizes = SAMTOOLS_FAIDX.out.sizes
    versions = versions.mix(SAMTOOLS_FAIDX.out.versions)

    /*
     * Convert BAM to BEDGraph
     */
    ch_sorted_bam
        .combine([1])
        .set { ch_genomecov_input }
    ch_sizes = ch_chr_sizes.map { meta, sizes -> sizes }
    extension = 'bedGraph'
    to_sort   = false
    BEDTOOLS_GENOMECOV ( ch_genomecov_input, ch_sizes, extension, to_sort )
    ch_bedgraph      = BEDTOOLS_GENOMECOV.out.genomecov
    versions = versions.mix(BEDTOOLS_GENOMECOV.out.versions)

    /*
     * Convert BEDGraph to BigWig
     */
    UCSC_BEDGRAPHTOBIGWIG ( ch_bedgraph, ch_sizes )
    ch_bigwig = UCSC_BEDGRAPHTOBIGWIG.out.bigwig
    versions = versions.mix(UCSC_BEDGRAPHTOBIGWIG.out.versions)

    emit:
    ch_bigwig
    ch_bedgraph
    versions
}
