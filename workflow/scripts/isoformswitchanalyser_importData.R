
source(".Rprofile")
suppressPackageStartupMessages({
    library("logger")
    library("readr")
    library("dplyr")
    library("IsoformSwitchAnalyzeR")
    source("R/readGtf.R")
})

BPPARAM = BiocParallel::MulticoreParam(snakemake@threads[[1]])

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

design_matrix <- sample_info %>%
    dplyr::select(sampleID = "sample_id",
        condition = condition_temp) %>%
    as.data.frame()


log_info("Import counts...")
if(snakemake@wildcards$dataset == "illumina") {

    samples <- snakemake@input$counts

    order <- match(sample_info$illumina,
        snakemake@input$counts %>%
            sub("/quant.sf$", "", .) %>%
            basename())

    names(samples) <- sample_info$sample_id[order]

    salmon_quant <- importIsoformExpression(
        sampleVector = samples)

    if(grepl("flair", snakemake@wildcards$annotation)) {
        # remove the gene name from the flair isoform ids
        salmon_quant$abundance$isoform_id <- salmon_quant$abundance$isoform_id %>%
            sub("_.+$", "", .)

        salmon_quant$counts$isoform_id <- salmon_quant$counts$isoform_id %>%
            sub("_.+$", "", .)
    }

} else if(snakemake@wildcards$dataset == "cdna"){

    counts <- snakemake@input$counts %>%
        purrr::set_names(basename(.) %>%
                             sub("_quant.tsv", "", .)) %>%
        purrr::map(read_tsv, col_types = "ci") %>%
        bind_rows(.id = "barcode") %>%
        tidyr::pivot_wider(names_from = "barcode", values_from = NumReads,
            values_fill = 0)

    if(grepl("flair", snakemake@wildcards$annotation)) {
        # remove the gene name from the flair isoform ids
        counts$Name <- counts$Name %>%
            sub("_.+$", "", .)
    }

    if(grepl("ref", snakemake@wildcards$annotation)) {
        order <- match(sample_info$cdna, sub("cdna_merged_", '', names(counts))) - 1
    } else {
        order <- match(sample_info$cdna, names(counts)) - 1
    }

    names(counts) <- c("isoform_id", sample_info$sample_id[order])


    missing_transcripts <- read_gtf(snakemake@input$gtf) %>%
        filter(feature == "transcript") %>%
        pull(attributes) %>%
        sub('^.+transcript_id \"', "", .) %>%
        sub('\";.*$', "", .)


    # add empty rows for the sake of isoformswitchanalyser
    missing_transcripts <- missing_transcripts[!missing_transcripts %in% counts$isoform_id]

    counts <- tibble(isoform_id = missing_transcripts,
            "190220_2_iBAT" = 0,
            "190220_4_iBAT" = 0,
            "190220_9_iBAT" = 0,
            "190220_11_iBAT" = 0,
            "190220_14_iBAT" = 0,
            "190220_15_iBAT" = 0)  %>%
        bind_rows(counts)
}


log_info("Create switchlist...")
if(snakemake@wildcards$dataset == "illumina") {
        switchList <- importRdata(
        isoformCountMatrix   = salmon_quant$counts,
        isoformRepExpression = salmon_quant$abundance,
        designMatrix         = design_matrix,
        isoformExonAnnoation = snakemake@input$gtf,
        isoformNtFasta       = snakemake@input$transcripts)

} else {

    # in the isoformswitchanalyser vignette the isoformRepExpression is empty
    # no idea, what it calculates here
    # does it normalise for library size?
    # does it normalise for transcript length? (that would be wrong for nanopore data)
    switchList <- importRdata(
        isoformCountMatrix   = counts,
        designMatrix         = design_matrix,
        isoformExonAnnoation = snakemake@input$gtf,
        isoformNtFasta       = snakemake@input$transcripts)

}


log_info("Prefiltering")
switchList <- preFilter(switchList,
    geneExpressionCutoff = 3,
    isoformExpressionCutoff = 1,
    IFcutoff = .025,
    removeSingleIsoformGenes = TRUE,
    reduceToSwitchingGenes = FALSE
)


log_info("Testing for DTU using DRIMSeq")
switchList <- isoformSwitchTestDRIMSeq(switchList,
        alpha = snakemake@params$dtu_cutoff,
        dIFcutoff = snakemake@params$dIF_cutoff,
        testIntegration = "intersect",
        reduceToSwitchingGenes = FALSE,
        dmFilterArgs = list(
            min_samps_feature_expr = 3,
            min_feature_expr = snakemake@params$min_feature_expr,
            min_samps_feature_prop = 3,
            min_feature_prop = snakemake@params$min_feature_prop,
            min_samps_gene_expr = 6,
            min_gene_expr = snakemake@params$min_gene_expr),
        dmPrecisionArgs = list(BPPARAM = BPPARAM),
        dmFitArgs = list(BPPARAM = BPPARAM),
        dmTestArgs = list(BPPARAM = BPPARAM))


log_info("Exporting DTU results table...")
# this table holds all values, significant or not
switchList$isoformFeatures %>%
    as_tibble() %>%
    dplyr::select(gene_id,
        transcript_id = isoform_id,
        mgi_symbol = gene_name,
        gene_biotype,
        gene = gene_switch_q_value,
        transcript = isoform_switch_q_value,
        proportion_22 = IF1,
        proportion_4 = IF2) %>%
    write_csv(snakemake@output$res)


if(snakemake@wildcards$annotation != "ref") {
    log_info("Annotating ORFs...")
    switchList <- analyzeORF(switchList)
}


log_info("Add aa sequence...")
switchList <- extractSequence(switchList,
    onlySwitchingGenes = TRUE,
    alpha = snakemake@params$dtu_cutoff,
    dIFcutoff = snakemake@params$dIF_cutoff,
    extractNTseq = FALSE,
    writeToFile = TRUE,
    pathToOutput = snakemake@params[["aa_path"]],
    outputPrefix = snakemake@params[["aa_prefix"]])


log_info("Predict alternative splicing and intron retention...")
switchList <- analyzeAlternativeSplicing(switchList,
    onlySwitchingGenes = TRUE,
    alpha = snakemake@params$dtu_cutoff,
    dIFcutoff = snakemake@params$dIF_cutoff)

switchList <- analyzeIntronRetention(switchList,
    onlySwitchingGenes = TRUE,
    alpha = snakemake@params$dtu_cutoff,
    dIFcutoff = snakemake@params$dIF_cutoff)


log_info("Analyze orf seq similarity...")
switchList <- analyzeSwitchConsequences(switchList,
    c("tss", "tts", "intron_retention", "ORF_seq_similarity", "NMD_status"),
    alpha = snakemake@params$dtu_cutoff,
    dIFcutoff = snakemake@params$dIF_cutoff)

log_info("Writing to disc...")
saveRDS(switchList, snakemake@output[[1]])

log_success("Done.")
