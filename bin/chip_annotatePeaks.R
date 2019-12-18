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
        bindingRegion = c(-2000, 1000)) %>%
    map(as_tibble) %>%
    map(select, peak, tx_name) %>%
    imap(.x = ., ~ set_names(.x, c(.y, "tx_name"))) %>%
    purrr::reduce(full_join, by = "tx_name")

# dataframe of every transcript and wether it has a peak or not at its TSS
transcripts %>%
    as_tibble() %>%
    left_join(peaksAtTSS, by = "tx_name") %>%
    mutate_at(vars(ends_with(".bed")),
        function(x) if_else(is.na(x), FALSE, TRUE)) %>%
    dplyr::rename(ensembl_transcript_id_version = "tx_name",
        ensembl_gene_id_version = "gene_id") %>%
    group_by(ensembl_transcript_id_version) %>%
    summarise_at(vars(matches(".bed")), any) %>%
    write_csv(snakemake@output[["TSS"]])

# annotate peaks overlapping a TSS, on gene level
peaksAtTSS <- gr %>%
    map(annoPeaks, annoData = genes,
        bindingType = "startSite",
        bindingRegion = c(-2000, 1000)) %>%
    map(as_tibble) %>%
    map(select, peak, feature) %>%
    imap(.x = ., ~ set_names(.x, c(paste0(.y, "_TSS"), "feature"))) %>%
    purrr::reduce(full_join, by = "feature")

# annotate peaks overlapping a gene
allPeaks <- gr %>%
    map(annoPeaks, annoData = genes,
        bindingType = "fullRange",
        bindingRegion = c(-2000, 1000)) %>%
    map(as_tibble) %>%
    map(select, peak, feature) %>%
    imap(.x = ., ~ set_names(.x, c(paste0(.y, "_fullRange"), "feature"))) %>%
    purrr::reduce(full_join, by = "feature")

# remove peaks overlapping a transcript
internal_NW <- allPeaks %>%
    filter(!NW_broad.bed_fullRange %in% peaksAtTSS$NW_broad.bed_TSS) %>%
    pull(NW_broad.bed_fullRange)

internal_NC <- allPeaks %>%
    filter(!NC_broad.bed_fullRange %in% peaksAtTSS$NC_broad.bed_TSS) %>%
    pull(NC_broad.bed_fullRange)


# peaks that are in full range but not at TSS should be internal
df <- genes@ranges %>%
    as_tibble() %>%
    dplyr::rename(feature = "names") %>%
    left_join(peaksAtTSS, by = "feature") %>%
    left_join(allPeaks, by = "feature") %>%
    mutate(
        NW_broad.bed_internal = if_else(NW_broad.bed_fullRange %in% internal_NW,
                                        NW_broad.bed_fullRange, NA_character_),
        NC_broad.bed_internal = if_else(NC_broad.bed_fullRange %in% internal_NC,
                                        NC_broad.bed_fullRange, NA_character_))

df1 <- df %>%
    dplyr::select(-NC_broad.bed_fullRange, -NW_broad.bed_fullRange) %>%
    group_by(feature) %>%
    summarise_at(vars(matches(".bed_")), function(x) any(!is.na(x))) %>%
    dplyr::rename(ensembl_gene_id_version = "feature")

df2 <- df %>%
    group_by(feature) %>%
    summarise_at(vars(matches(".bed_fullRange")),
        function(x) length(unique(na.omit(x)))) %>%
    dplyr::rename(ensembl_gene_id_version = "feature")


df1 %>%
    left_join(df2, by = "ensembl_gene_id_version") %>%
    write_csv(snakemake@output[["full_range"]])
