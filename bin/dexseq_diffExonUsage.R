#!/usr/bin/Rscript --no-restore --no-environ --no-save

# set libpaths to packrat local library
source("packrat/init.R")
suppressPackageStartupMessages(library("DEXSeq"))
    select <- dplyr::select
BPPARAM = BiocParallel::MulticoreParam(snakemake@threads[[1]])

################################################################################
#
# author: Christoph Kiefer
# email: christophak@bmb.sdu.dk
#
################################################################################

message("Importing data...")
dxd <- readRDS(snakemake@input[[1]])

message("Normalising...")
dxd <- estimateSizeFactors(dxd)

message("Calculating dispersion estimates...")
dxd <- estimateDispersions(dxd, BPPARAM = BPPARAM)

message("Testing for differential exon usage...")
dxd <- testForDEU(dxd, BPPARAM = BPPARAM)
dxd <- estimateExonFoldChanges(dxd, fitExpToVar = "condition_temp",
    BPPARAM = BPPARAM)

message("Saving results to disc...")
saveRDS(dxd, snakemake@output[["dxd"]])

message("Generating report...")
dxr <- DEXSeqResults(dxd)

ensembl <- biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")

DEXSeqHTML(dxr,
    path = dirname(snakemake@output[["report"]]),
    file = basename(snakemake@output[["report"]]),
    fitExpToVar = "condition_temp",
    FDR = .05,
    mart = ensembl,
    filter = "ensembl_gene_id_version",
    attributes = c("mgi_symbol","description"),
    BPPARAM = BPPARAM)

message("Done")
