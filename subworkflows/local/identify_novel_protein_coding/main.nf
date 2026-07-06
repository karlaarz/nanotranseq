include { CPAT_BUILD_MODEL          } from '../../../modules/local/cpat_build_model'
include { CPAT                      } from '../../../modules/local/cpat'
include { PLEK_RUN                  } from '../../../modules/local/plek_run'
include { PLEK_PARSE                } from '../../../modules/local/plek_parse'
include { FEELNC_CODPOT_RUN         } from '../../../modules/local/feelnc_codpot_run'
include { FEELNC_CODPOT_PARSE       } from '../../../modules/local/feelnc_codpot_parse'
include { COMBINE_PREDICTIONS       } from '../../../modules/local/combine_predictions'


workflow IDENTIFY_NOVEL_PROTEIN_CODING {
    take:
        ch_novel_tmap                         // NOVEL_TRANSCRIPTS.out.tmap
        ch_novel_gtf                          // NOVEL_TRANSCRIPTS.out.novel_gtf
        ch_novel_fasta                        // NOVEL_TRANSCRIPTS.out.novel_fasta
        ch_cds_fasta                          // coding training FASTA for CPAT
        ch_noncoding_fasta                    // noncoding reference FASTA for CPAT and FEELnc
        ch_mrna_fasta                         // protein-coding mRNA reference FASTA for FEELnc


    main:
        ch_versions = Channel.empty()

        ch_candidate_gtf = ch_novel_gtf
        ch_candidate_fa  = ch_novel_fasta
        ch_tmap          = ch_novel_tmap.map { meta, tmap -> tmap }

        //
        // CPAT: Prepare or build models
        //
        ch_cpat_hexamer = Channel.empty()
        ch_cpat_logit = Channel.empty()

        if (!params.skip_cpat && params.cpat_hexamer && params.cpat_logit_model) {
            // Use provided models
            ch_cpat_hexamer = Channel.fromPath(params.cpat_hexamer)
            ch_cpat_logit = Channel.fromPath(params.cpat_logit_model)
        } else if (!params.skip_cpat) {
            // Build models
            CPAT_BUILD_MODEL (
                ch_cds_fasta, //coding sequences only
                ch_noncoding_fasta
            )
            ch_cpat_hexamer = CPAT_BUILD_MODEL.out.hexamer
            ch_cpat_logit = CPAT_BUILD_MODEL.out.logit_model
            ch_versions = ch_versions.mix(CPAT_BUILD_MODEL.out.versions)
        }

        //
        // CPAT: Coding Potential Assessment Tool
        //
        if (!params.skip_cpat) {
            CPAT (
                ch_candidate_fa,
                ch_cpat_hexamer,
                ch_cpat_logit
            )
            ch_cpat_results = CPAT.out.cpat_results
            ch_versions = ch_versions.mix(CPAT.out.versions)
        }

        //
        // FEELnc: FlExible Extraction of LncRNAs
        //
        if (!params.skip_feelnc) {
            FEELNC_CODPOT_RUN (
                ch_candidate_fa,
                ch_mrna_fasta,
                ch_noncoding_fasta
            )

            FEELNC_CODPOT_PARSE (
                FEELNC_CODPOT_RUN.out.codpot_full,
                ch_candidate_fa
            )

            ch_feelnc_results = FEELNC_CODPOT_PARSE.out.feelnc_results
            ch_versions = ch_versions.mix(FEELNC_CODPOT_RUN.out.versions)
            ch_versions = ch_versions.mix(FEELNC_CODPOT_PARSE.out.versions)
        }

        //
        // PLEK: Predictor of lncRNAs and mRNAs based on k-mer
        //
        if (!params.skip_plek) {
            PLEK_RUN (
                ch_candidate_fa
            )

            PLEK_PARSE (
                PLEK_RUN.out.plek_raw
            )
            ch_plek_results = PLEK_PARSE.out.plek_results
            ch_versions = ch_versions.mix(PLEK_RUN.out.versions)
            ch_versions = ch_versions.mix(PLEK_PARSE.out.versions)
        }

        //
        // COMBINE: Merge predictions from CPAT, FEELnc and PLEK
        //
        COMBINE_PREDICTIONS (
            ch_cpat_results,
            ch_feelnc_results,
            ch_plek_results,
            ch_candidate_gtf,
            ch_candidate_fa,
            ch_tmap
        )
        ch_final_protein_coding_gtf = COMBINE_PREDICTIONS.out.protein_coding_gtf
        ch_final_protein_coding_fa  = COMBINE_PREDICTIONS.out.protein_coding_fasta
        ch_protein_coding_pred_summary = COMBINE_PREDICTIONS.out.protein_coding_pred_summary
        ch_coding_potential_report = COMBINE_PREDICTIONS.out.report
        ch_versions = ch_versions.mix(COMBINE_PREDICTIONS.out.versions)


    emit:
        candidate_gtf                   = ch_candidate_gtf                   // Novel transcript candidate GTF
        candidate_fasta                 = ch_candidate_fa                    // Novel transcript candidate FASTA
        cpat_hexamer                    = ch_cpat_hexamer                    // CPAT hexamer (provided or built)
        cpat_logit                      = ch_cpat_logit                      // CPAT logit (provided or built)
        cpat_results                    = ch_cpat_results                    // CPAT predictions
        feelnc_results                  = ch_feelnc_results                  // FEELnc predictions
        plek_results                    = ch_plek_results                    // PLEK predictions
        final_protein_coding_gtf        = ch_final_protein_coding_gtf        // Final protein-coding GTF
        final_protein_coding_fasta      = ch_final_protein_coding_fa         // Final protein-coding FASTA
        report                          = ch_coding_potential_report         // Summary report
        protein_coding_pred_summary     = ch_protein_coding_pred_summary     // protein-coding prediction
        versions                        = ch_versions
}