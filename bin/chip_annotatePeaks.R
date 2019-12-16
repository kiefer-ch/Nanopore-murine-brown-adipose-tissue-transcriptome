#!/usr/bin/Rscript --no-restore --no-environ --no-save

source("packrat/init.R")
library("ChIPpeakAnno")
library("AnnotationDbi")
library("dplyr")
library("purrr")
library("readr")

################################################################################
#
# author: Christoph Kiefer
# email: christophak@bmb.sdu.dk
#
################################################################################

# import txdb
txdb <- loadDb(snakemake@input[["txdb"]])
genes <- toGRanges(txdb, feature = "gene")
transcripts <- toGRanges(txdb, feature = "transcript")

# import macs2 bedfiles
gr <- c(snakemake@input[["h3k4_warm"]], snakemake@input[["h3k4_cold"]]) %>%
    set_names(basename(.)) %>%
    map(toGRanges, format = "BED", header = FALSE)

# annotate peaks overlapping a TSS, on tx level
peaksAtTSS <- gr %>%
    map(annoPeaks, annoData = transcripts,
        bindingType = "startSite",
        bindingRegion = c(-2000, 500)) %>%
    map(as_tibble) %>%
    map(select, peak, tx_name) %>%
    imap(.x = ., ~ set_names(.x, c(.y, "tx_name"))) %>%
    purrr::reduce(left_join, by = "tx_name")

# dataframe of every transcript and wether it has a peak or not at its TSS
transcripts %>%
    as_tibble() %>%
    left_join(peaksAtTSS, by = "tx_name") %>%
    mutate_at(vars(ends_with(".bed")),
        function(x) if_else(is.na(x), FALSE, TRUE)) %>%
    dplyr::rename(ensembl_transcript_id_version = "tx_name",
        ensembl_gene_id_version = "gene_id") %>%
    write_csv(snakemake@output[["TSS"]])

# annotate peaks overlapping a TSS, on gene level
peaksAtTSS <- gr %>%
    map(annoPeaks, annoData = genes,
        bindingType = "startSite",
        bindingRegion = c(-2000, 500)) %>%
    map(as_tibble) %>%
    map(select, peak, feature) %>%
    imap(.x = ., ~ set_names(.x, c(paste0(.y, "_TSS"), "feature"))) %>%
    purrr::reduce(left_join, by = "feature")

# annotate peaks overlapping a gene
allPeaks <- gr %>%
    map(annoPeaks, annoData = genes,
        bindingType = "fullRange",
        bindingRegion = c(-2000, 500)) %>%
    map(as_tibble) %>%
    map(select, peak, feature) %>%
    imap(.x = ., ~ set_names(.x, c(paste0(.y, "_fullRange"), "feature"))) %>%
    purrr::reduce(left_join, by = "feature")

# peaks that are in full range but not at TSS should be internal
genes@ranges %>%
    as_tibble() %>%
    dplyr::rename(feature = "names") %>%
    left_join(peaksAtTSS) %>%
    left_join(allPeaks) %>%
    dplyr::rename(ensembl_gene_id_version = "feature") %>%
    mutate_at(vars(matches(".bed")),
        function(x) if_else(is.na(x), FALSE, TRUE)) %>%
    write_csv(snakemake@output[["full_range"]])
