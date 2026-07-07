//
// RAW READS QUALITY CONTROL WORKFLOW
//

include { FASTQC      } from '../../../modules/nf-core/fastqc/main'
include { TOULLIGQC   } from '../../../modules/nf-core/toulligqc/main'
include { NANOPLOT    } from '../../../modules/nf-core/nanoplot/main'
include { MULTIQC     } from '../../../modules/nf-core/multiqc/main'

workflow RAW_READS_QC {
    take:
    reads    // Raw reads input channel

    main:
    versions = Channel.empty()

    // Run FASTQC
    FASTQC(reads)
    fastqc_zip = FASTQC.out.zip
    fastqc_html = FASTQC.out.html

    // MULTIQC(FASTQC.out.zip.collect{it[1]},
    //         [],
    //         [],
    //         [],
    //         [],
    //         [])

    MULTIQC(
          FASTQC.out.zip
              .collect { it[1] }
              .map { files -> tuple([:], files, [], [], [], []) }
    )

    // Run TOULLIGQC
    TOULLIGQC(reads)
    toulligqc_report_data = TOULLIGQC.out.report_data
    toulligqc_report_html = TOULLIGQC.out.report_html
    toulligqc_plots_html  = TOULLIGQC.out.plots_html
    toulligqc_plotly_js   = TOULLIGQC.out.plotly_js

    // Run Nanoplot
    NANOPLOT(reads)
    nanoplot_html = NANOPLOT.out.html
    nanoplot_png  = NANOPLOT.out.png
    nanoplot_txt  = NANOPLOT.out.txt

    // Collect versions for all tools used in this workflow
    //versions = versions.mix(FASTQC.out.versions, MULTIQC.out.versions, TOULLIGQC.out.versions, NANOPLOT.out.versions)

    emit:
    fastqc_zip
    fastqc_html

    multiqc_report = MULTIQC.out.report

    toulligqc_report_data
    toulligqc_report_html
    toulligqc_plots_html
    toulligqc_plotly_js

    nanoplot_html
    nanoplot_png
    nanoplot_txt

    versions
}
