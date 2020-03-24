#!/usr/bin/Rscript --no-restore --no-environ --no-save

source(".Rprofile")
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("DESeq2"))

################################################################################
#
# author: Christoph Kiefer
# email: christophak@bmb.sdu.dk
#
################################################################################

# tx and gene level
if (snakemake@params[["level"]] == "txlevel") {
    id <- "ensembl_transcript_id_version"
} else if (snakemake@params[["level"]] == "genelevel") {
    id <- "ensembl_gene_id_version"
}

# import data
dds <- readRDS(snakemake@input[["dds"]])
biomart <- readRDS(snakemake@input[["biomart"]])

# prepare countmatrices
if (snakemake@params[["vst"]] == 1) {
    rld <- vst(dds, blind = FALSE)
} else {
    rld <- rlog(dds, blind = FALSE)
}

ntd <- normTransform(dds)


# with variance shrinking
assay(rld) %>%
    as_tibble(rownames = id) %>%
    left_join(biomart, by = id) %>%
    dplyr::select(starts_with("ensembl_"), mgi_symbol, description,
        gene_biotype, everything()) %>%
    write_csv(snakemake@output[["rld"]])

# without
assay(ntd) %>%
    as_tibble(rownames = id) %>%
    left_join(biomart, by = id) %>%
    dplyr::select(starts_with("ensembl_"), mgi_symbol, description,
        gene_biotype, everything()) %>%
    write_csv(snakemake@output[["ntd"]])

# tpm
# https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/
fpkmToTpm <- function(fpkm) {
    exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

if (as.logical(snakemake@params[["tpm"]])) {
    tpm <- fpkm(dds, robust = TRUE)
    tpm <- apply(tpm, 2, fpkmToTpm)

    tpm %>%
        as_tibble(rownames = id) %>%
        left_join(biomart, by = id) %>%
        dplyr::select(starts_with("ensembl_"), mgi_symbol, description,
            gene_biotype, everything()) %>%
        write_csv(snakemake@output[["tpm"]])
}

# raw counts
counts(dds, normalized = FALSE) %>%
    as_tibble(rownames = id) %>%
    left_join(biomart, by = id) %>%
    dplyr::select(starts_with("ensembl_"), mgi_symbol, description,
        gene_biotype, everything()) %>%
    write_csv(snakemake@output[["cts"]])
