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

# get gene ids and annotation
if(file.exists(snakemake@output[["gene"]])) {
    message("Skipping gene annotation.")
} else {
message("Fetching gene annotation...")
gene_ids <- keys(txdb, "GENEID")

annotation <- biomaRt::getBM(
        attributes = c("ensembl_gene_id_version",
            "mgi_symbol", "description",
            "gene_biotype"),
        filters = "ensembl_gene_id_version",
        values = unique(gene_ids),
        mart = ensembl) %>%
    as_tibble() %>%
    distinct(ensembl_gene_id_version, .keep_all = TRUE) %>%
    saveRDS(snakemake@output[["gene"]])
}
