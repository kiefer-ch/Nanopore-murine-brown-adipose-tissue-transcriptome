
source(".Rprofile")
suppressPackageStartupMessages({
    library("logger")
    library("readr")
    library("dplyr")
    library("IsoformSwitchAnalyzeR")
    source("R/readGtf.R")
})

################################################################################
#
# author: Christoph Kiefer
# email: christophak@bmb.sdu.dk
#
################################################################################

# prepare sample_info
sample_info <- read_csv(snakemake@input[["sample_info"]],
        show_col_types = FALSE) %>%
    filter(!is.na(ont)) %>%
    mutate_at(vars(matches("condition")), as.factor)

# get design matrix
if (snakemake@wildcards$annotation %in% c("rna_flair", "rna_stringtie")) {
    design_matrix <- data.frame(
            sampleID = c("rt", "cool"),
            condition = as.factor(c("22", "4")))
} else {
    design_matrix <- sample_info %>%
        dplyr::select(
            sampleID = "sample_id",
            condition = condition_temp) %>%
        as.data.frame()
}

# get count matrix
all_transcripts <- read_gtf(snakemake@input$gtf) %>%
    filter(feature == "transcript") %>%
    pull(attributes) %>%
    sub('^.+transcript_id \"', "", .) %>%
    sub('\";.*$', "", .)

if (snakemake@wildcards$annotation %in% c("rna_flair", "rna_stringtie")) {
    counts = tibble(
        isoform_id = all_transcripts,
        rt = 100,
        cool = 100)
} else {
    counts = tibble(
        isoform_id = all_transcripts,
        "190220_2_iBAT" = runif(length(all_transcripts), 100, 1000),
        "190220_4_iBAT" = runif(length(all_transcripts), 100, 1000),
        "190220_9_iBAT" = runif(length(all_transcripts), 100, 1000),
        "190220_11_iBAT" = runif(length(all_transcripts), 100, 1000),
        "190220_14_iBAT" = runif(length(all_transcripts), 100, 1000),
        "190220_15_iBAT" = runif(length(all_transcripts), 100, 1000))
}


# make switchlist
switchList <- importRdata(
    isoformCountMatrix   = counts,
    designMatrix         = design_matrix,
    isoformExonAnnoation = snakemake@input$gtf,
    isoformNtFasta = snakemake@input$transcripts)


log_info("Predict alternative splicing and intron retention...")
switchList <- analyzeAlternativeSplicing(switchList,
    onlySwitchingGenes = FALSE)

switchList <- analyzeIntronRetention(switchList,
    onlySwitchingGenes = FALSE,
    alpha = 1,
    dIFcutoff = 0)

log_info("Writing to disc...")
saveRDS(switchList, snakemake@output[[1]])

log_success("Done.")
