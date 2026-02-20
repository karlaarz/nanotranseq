#!/usr/bin/env Rscript

library(IsoformSwitchAnalyzeR)

# Nextflow template variables
design_file <- "${design_file}"
fasta_file <- "${fasta}"
gtf_file <- "${gtf}"
gene_results_file <- "${gene_results_file}"
transcript_results_file <- "${transcript_results_file}"

# Read design file
samples <- read.csv(design_file, stringsAsFactors=FALSE, header=TRUE)
colnames(samples) <- c("sampleID", "condition")

# Find count and TPM files
count_files <- list.files(pattern = "transcript_counts.tsv")
tpm_files <- list.files(pattern = "transcript_tpm.tsv")

# Helper to read data (handles both merged and individual files)
read_data <- function(files, samples, type="counts") {
    # Check for merged file
    merged_pattern <- paste0("all_samples.*", type, ".tsv")
    merged_file <- grep(merged_pattern, files, value=TRUE)

    if (length(merged_file) > 0) {
        message("Detected merged file for ", type, ": ", merged_file[1])
        d <- read.table(merged_file[1], header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names=FALSE)

        # Ensure first two columns are isoform_id and gene_id
        colnames(d)[1] <- "isoform_id"
        colnames(d)[2] <- "gene_id"

        # Check if all samples are present
        missing_samples <- setdiff(samples\$sampleID, colnames(d))
        if (length(missing_samples) > 0) {
            stop("Merged file ", merged_file[1], " is missing columns for samples: ", paste(missing_samples, collapse=", "))
        }

        # Select and order columns: isoform_id, sample1, sample2...
        # importRdata expects only isoform_id and counts (no gene_id column)
        d <- d[, c("isoform_id", samples\$sampleID)]

        # Ensure IDs are character
        d\$isoform_id <- as.character(d\$isoform_id)

        # Force count columns to be numeric
        for (col in samples\$sampleID) {
            d[[col]] <- as.numeric(as.character(d[[col]]))
        }

        message("Debug: Column classes for ", type)
        # print(sapply(d, class))

        return(d)
    }

    # Fallback to individual files
    matched_files <- sapply(samples\$sampleID, function(sid) {
        # Use [.] instead of \\. to avoid escape issues in Nextflow templates
        sid_regex <- paste0("^", sid, ".*", type, "[.]tsv\$")
        match <- grep(sid_regex, files, value= TRUE)
        if (length(match) == 0) {
             match <- grep(sid, files, value= TRUE)
        }
        if (length(match) == 0) return(NA)
        return(match[1])
    })

    if (any(is.na(matched_files))) {
        stop("Missing files for ", type, ": ", paste(samples\$sampleID[is.na(matched_files)], collapse=", "))
    }

    first <- read.table(matched_files[1], header=TRUE, sep="\t", stringsAsFactors=FALSE)
    # Col 1: tx, Col 2: gene (ignore gene_id for matrix)
    df <- data.frame(isoform_id = first[, 1])

    for (i in 1:nrow(samples)) {
        sid <- samples\$sampleID[i]
        f <- matched_files[i]
        d <- read.table(f, header=TRUE, sep="\t", stringsAsFactors=FALSE)
        if (!identical(d[,1], df\$isoform_id)) {
            d <- d[match(df\$isoform_id, d[,1]), ]
        }
        df[[sid]] <- as.numeric(d[, 3])
    }
    return(df)
}

tpm_df <- read_data(tpm_files, samples, "transcript_tpm")
counts_df <- read_data(count_files, samples, "transcript_counts")

# Create switchAnalyzeRlist
switchList <- importRdata(
    isoformCountMatrix = counts_df,
    isoformRepExpression = tpm_df,
    designMatrix = samples,
    isoformExonAnnoation = gtf_file,
    showProgress = FALSE,
    ignoreAfterPeriod = TRUE
)

# Add sequences manually using Biostrings to avoid dependency issues with addIsoformNtSequence
message("Adding sequences using Biostrings...")

if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Biostrings package is required but not installed.")
}

dna <- Biostrings::readDNAStringSet(fasta_file)
# Clean IDs (remove description after space)
names(dna) <- sub("[[:space:]].*", "", names(dna))

# Match sequences with isoform IDs
    common <- intersect(names(dna), switchList\$isoformFeatures\$isoform_id)

    if(length(common) == 0) {
        # Try removing versions (text after last dot)
        names(dna) <- sub("[.][^.]*\$", "", names(dna))
        common <- intersect(names(dna), switchList\$isoformFeatures\$isoform_id)
    }

    if(length(common) == 0) {
         # Try removing pipes (text after first pipe)
         names(dna) <- sub("[|].*", "", names(dna))
         common <- intersect(names(dna), switchList\$isoformFeatures\$isoform_id)
    }

if(length(common) > 0) {
    # Assign sequences to the switchList object
    switchList\$ntSequence <- dna[common]
    message(paste("Successfully added", length(common), "sequences."))
} else {
    warning("No matching sequences found between FASTA and isoform IDs.")
}

# Run analysis pipeline with error handling
tryCatch({
    switchList <- preFilter(switchList)

    # Import external results (DRIMSeq or DEXSeq)
    # We expect these files to be present if we are running in this mode

    message("Importing external statistical test results...")

    # Check if files exist and are not empty
    if (file.exists(gene_results_file) && file.exists(transcript_results_file)) {

        # Determine if it is DRIMSeq or DEXSeq based on column names or file name logic?
        # Actually, we can just look for the expected columns: gene_id, feature_id, adj_pvalue

        gene_res <- read.csv(gene_results_file, stringsAsFactors=FALSE)
        tx_res <- read.csv(transcript_results_file, stringsAsFactors=FALSE)

        if (nrow(gene_res) > 0 && nrow(tx_res) > 0) {

             # Map gene q-values
             # Standardize columns
             if (!"gene_id" %in% colnames(gene_res) || !"adj_pvalue" %in% colnames(gene_res)) {
                 warning("Gene results file missing expected columns (gene_id, adj_pvalue)")
             } else {
                 gene_q <- gene_res[, c("gene_id", "adj_pvalue")]
                 colnames(gene_q) <- c("gene_id", "gene_switch_q_value")
                 switchList\$isoformFeatures <- merge(switchList\$isoformFeatures, gene_q, by="gene_id", all.x=TRUE)
             }

             # Map transcript q-values
             if (!"feature_id" %in% colnames(tx_res) || !"adj_pvalue" %in% colnames(tx_res)) {
                 warning("Transcript results file missing expected columns (feature_id, adj_pvalue)")
             } else {
                 tx_q <- tx_res[, c("feature_id", "adj_pvalue")]
                 colnames(tx_q) <- c("isoform_id", "isoform_switch_q_value")
                 switchList\$isoformFeatures <- merge(switchList\$isoformFeatures, tx_q, by="isoform_id", all.x=TRUE)
             }

             # Fill NAs
             switchList\$isoformFeatures\$isoform_switch_q_value[is.na(switchList\$isoformFeatures\$isoform_switch_q_value)] <- 1
             switchList\$isoformFeatures\$gene_switch_q_value[is.na(switchList\$isoformFeatures\$gene_switch_q_value)] <- 1

        } else {
            message("External results are empty. Skipping switch identification.")
            switchList\$isoformFeatures\$isoform_switch_q_value <- 1
            switchList\$isoformFeatures\$gene_switch_q_value <- 1
        }
    } else {
        message("Result files not found. Skipping switch identification.")
        switchList\$isoformFeatures\$isoform_switch_q_value <- 1
        switchList\$isoformFeatures\$gene_switch_q_value <- 1
    }

    # Analysis
    switchList <- analyzeORF(switchList, orfMethod = "longest")
    switchList <- extractSequence(switchList)
    switchList <- analyzeAlternativeSplicing(switchList)

    # Summary and Plots
    switchList <- switchPlot(switchList)

    saveRDS(switchList, "isoform_switch_results.rds")

}, error = function(e) {
    message("IsoformSwitchAnalyzeR analysis failed.")
    message(e)
    # Create empty output to prevent pipeline crash
    write.csv(data.frame(), "isoform_switch_results.csv")
    pdf("isoform_switch_plots.pdf")
    plot(1, 1, main="Analysis Failed")
    dev.off()
    saveRDS(switchList, "isoform_switch_results.rds")
})

# Versions
r_version <- paste(R.version\$major, R.version\$minor, sep=".")
isoform_version <- as.character(packageVersion('IsoformSwitchAnalyzeR'))

versions_text <- paste0(
    '"', "${task.process}", '":\n',
    '    r-base: ', r_version, '\n',
    '    bioconductor-isoformswitchanalyzer: ', isoform_version, '\n'
)
writeLines(versions_text, "versions.yml")
