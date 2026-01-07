include { PARABRICKS_MINIMAP2 } from '../../../modules/nf-core/parabricks/minimap2/main'
include { MINIMAP2_INDEX      } from '../../../modules/nf-core/minimap2/index/main'
include { MINIMAP2_ALIGN      } from '../../../modules/nf-core/minimap2/align/main'

workflow ALIGNMENT {

    take:
    use_gpus                // boolean: use gpus during alignment. default: false
    reads                   // reads channel: val(meta1),  path(reads)
    genome_reference        // genome refence: val(meta2),  path(fasta)
    minimap2_index          // minimap2 pre-built reference: path(reference)

    main:
    versions = Channel.empty()

    minimap2_bam = Channel.empty()
    minimap2_bai = Channel.empty()

    // Run Parabricks implementation of Minimap2 if use_gpus = true. 
    //If not, run standard implementation
    if (use_gpus == true) {

        interval = tuple([], [])
        known_sites = tuple([], [])
        output_fmt = 'bai'

        PARABRICKS_MINIMAP2(reads,
                            genome_reference,
                            interval,
                            known_sites,
                            output_fmt)

        minimap2_bam = PARABRICKS_MINIMAP2.out.bam
        minimap2_bai = PARABRICKS_MINIMAP2.out.bai

        versions = versions.mix(PARABRICKS_MINIMAP2.out.versions)

    } else {

        // If a Minimap2 index is provided, set fasta to null to skip indexing step
        index_input = minimap2_index
            .filter { it.name == 'no_minimap2_index' }
            .combine(genome_reference)
            .map { _, meta, fasta -> tuple(meta, fasta) }

        MINIMAP2_INDEX(index_input)

        final_index = MINIMAP2_INDEX.out.index
            .collect()

        // Align samples based on reference genome
        bam_format = true
        bam_index_extension = 'bai'
        cigar_paf_format = false
        cigar_bam = true

        MINIMAP2_ALIGN(reads,
                       final_index,
                       bam_format,
                       bam_index_extension,
                       cigar_paf_format,
                       cigar_bam)

        minimap2_bam = MINIMAP2_ALIGN.out.bam
        minimap2_bai = MINIMAP2_ALIGN.out.index

        versions = versions.mix(MINIMAP2_ALIGN.out.versions)

    }

    emit:
    versions

    minimap2_bam
    minimap2_bai
}
