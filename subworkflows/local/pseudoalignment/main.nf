include { SALMON_INDEX      } from '../../../modules/nf-core/salmon/index/main'
include { SALMON_QUANT      } from '../../../modules/nf-core/salmon/quant/main'
include { CUSTOM_TX2GENE    } from '../../../modules/nf-core/custom/tx2gene'
include { TXIMETA_TXIMPORT  } from '../../../modules/nf-core/tximeta/tximport'

workflow PSEUDOALIGNMENT {

    take:
    reads                           // channel: reads: val(meta), path(fastq)
    fasta                           // channel: fasta: val(meta), path(fasta)
    transcript_fasta                // channel: transcript_fasta: path(fasta)
    reference_gtf                   // channel: reference GTF: path(gtf)
    gene_id                         // channel: attributed gene ID: val
    gene_attributes                 // channel: extra gene attributes: val

    main:
    versions = Channel.empty()

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

    // Build tx2gene GTF reference
    CUSTOM_TX2GENE (
        reference_gtf.map { [ [id:"reference gtf"], it ] },
        salmon_out.collect{ it[1] }.map{ [ [id:"all_samples"], it] },
        Channel.value('salmon'),
        gene_id,
        gene_attributes,
    )
    versions = versions.mix(CUSTOM_TX2GENE.out.versions)

    // Parse quantification files from all samples
    TXIMETA_TXIMPORT (
        salmon_out
            .map{ it[1] }
            .collect()
            .map{ [ [id:"all_samples"], it] },
        CUSTOM_TX2GENE.out.tx2gene,
        Channel.value('salmon'),
    )
    versions = versions.mix(TXIMETA_TXIMPORT.out.versions)

    tpm_gene = TXIMETA_TXIMPORT.out.tpm_gene
    counts_gene = TXIMETA_TXIMPORT.out.counts_gene
    counts_gene_length_scaled = TXIMETA_TXIMPORT.out.counts_gene_length_scaled
    counts_gene_scaled = TXIMETA_TXIMPORT.out.counts_gene_scaled
    lengths_gene = TXIMETA_TXIMPORT.out.lengths_gene
    tpm_transcript = TXIMETA_TXIMPORT.out.tpm_transcript
    counts_transcript = TXIMETA_TXIMPORT.out.counts_transcript
    lengths_transcript = TXIMETA_TXIMPORT.out.lengths_transcript

    emit:
    versions

    tpm_gene
    counts_gene
    counts_gene_length_scaled
    counts_gene_scaled
    lengths_gene
    tpm_transcript
    counts_transcript
    lengths_transcript

}
