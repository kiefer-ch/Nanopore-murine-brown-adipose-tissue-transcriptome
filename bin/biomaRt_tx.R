#!/usr/bin/Rscript --no-restore --no-environ --no-save

# set libpaths to packrat local library
message("Loading packages...")
source("packrat/init.R")
suppressPackageStartupMessages(library("biomaRt"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("biomaRt"))
suppressPackageStartupMessages(library("AnnotationDbi"))

################################################################################
#
# author: Christoph Kiefer
# email: christophak@bmb.sdu.dk
#
################################################################################

ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
txdb <- loadDb(snakemake@input[["txdb"]])

# get transcript ids and annotations
message("Fetching transcript annotation...")

tx_ids <- keys(txdb, "TXNAME")

annotation <- biomaRt::getBM(
        attributes = c("ensembl_transcript_id_version",
            "ensembl_gene_id_version", "mgi_symbol", "description",
            "gene_biotype",
            "transcript_biotype", "transcript_length"),
        filters = "ensembl_transcript_id_version",
        values = unique(tx_ids),
        mart = ensembl) %>%
    as_tibble() %>%
    distinct(ensembl_transcript_id_version, .keep_all = TRUE) %>%
    saveRDS(snakemake@output[["tx"]])
