include { SAMTOOLS_FAIDX                } from '../../../modules/nf-core/samtools/faidx/main'
include { BEDTOOLS_GENOMECOV            } from '../../../modules/nf-core/bedtools/genomecov/main'
include { BEDTOOLS_BAMTOBED             } from '../../../modules/nf-core/bedtools/bamtobed/main'
include { UCSC_BEDGRAPHTOBIGWIG         } from '../../../modules/nf-core/ucsc/bedgraphtobigwig/main'
include { UCSC_BEDTOBIGBED              } from '../../../modules/nf-core/ucsc/bedtobigbed/main'

workflow BEDTOOLS_BIGWIG {
    take:
    ch_fasta
    ch_sorted_bam

    main:

    /*
     * Get chromosome sizes
     */
    SAMTOOLS_FAIDX(
        ch_fasta.map { meta, fasta -> tuple(meta, fasta, []) },
        true
    )
    ch_chr_sizes = SAMTOOLS_FAIDX.out.sizes

    /*
     * Convert BAM to BED
     */
    BEDTOOLS_BAMTOBED ( ch_sorted_bam )
    ch_bed = BEDTOOLS_BAMTOBED.out.bed

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

    /*
     * Convert BED to BigBed
     */
    UCSC_BEDTOBIGBED ( ch_bed, ch_sizes, [] )
    ch_bigbed = UCSC_BEDTOBIGBED.out.bigbed

    /*
     * Convert BEDGraph to BigWig
     */
    UCSC_BEDGRAPHTOBIGWIG ( ch_bedgraph, ch_sizes )
    ch_bigwig = UCSC_BEDGRAPHTOBIGWIG.out.bigwig

    emit:
    ch_bigwig
    ch_bedgraph
}
