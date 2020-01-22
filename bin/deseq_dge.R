#!/usr/bin/Rscript --no-restore --no-environ --no-save

source("packrat/init.R")
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("DESeq2"))
BPPARAM = BiocParallel::MulticoreParam(snakemake@threads[[1]])

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

# set reference
colData(dds)$condition_temp <- relevel(colData(dds)$condition_temp, ref = "22")

# filtering lowly expressed genes and fitting models
dds <- dds[rowSums(counts(dds)) > 10, ] %>%
    DESeq(., parallel = TRUE, BPPARAM = BPPARAM)

# calculating results table
res <- results(dds,
    contrast = c("condition_temp", "4", "22"),
    lfcThreshold = snakemake@params[["lfcThreshold"]],
    alpha = snakemake@params[["alpha"]],
    parallel = TRUE, BPPARAM = BPPARAM)

# shrinking log2 fold change
if (as.logical(snakemake@params[["shrink"]])) {
    res <- lfcShrink(dds,
        coef = "condition_temp_4_vs_22",
        res = res,
        lfcThreshold = snakemake@params[["lfcThreshold"]],
        type = "apeglm",
        parallel = TRUE, BPPARAM = BPPARAM)
}

# export
res %>%
    as_tibble(rownames = id) %>%
    left_join(biomart, by = id) %>%
    dplyr::select(starts_with("ensembl_"), mgi_symbol, description,
        gene_biotype, everything()) %>%
    write_csv(snakemake@output[[1]])
