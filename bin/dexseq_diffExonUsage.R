#!/usr/bin/Rscript --no-restore --no-environ --no-save

# set libpaths to packrat local library
source(".Rprofile")
suppressPackageStartupMessages(library("dplyr"))
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

message("Filtering lowly expressed genes...")
df <- featureCounts(dxd) %>%
    data.matrix() %>%
    rowMeans() %>%
    tibble::enframe() %>%
    tidyr::separate(name, c("gene", "exon"), sep = ':', remove = FALSE) %>%
    group_by(gene) %>%
    mutate(is_expressed = value > 5) %>%
    summarise(expressed_exons = sum(is_expressed))
n_tot <- nrow(df)
keep <- df %>%
    filter(expressed_exons >= 2) %>%
    pull(gene)
message(sprintf("Removed %s out of %s genes because of low expression.",
    n_tot - length(keep), n_tot))

keep <- rownames(dxd) %>%
    tibble::enframe() %>%
    tidyr::separate(value, c("gene", "exon"),
        sep = ':', remove = FALSE) %>%
    mutate(expressed = if_else(gene %in% keep, TRUE, FALSE)) %>%
    pull(expressed)

dxd <- dxd[keep, ]

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
    attributes = c("mgi_symbol", "description"),
    BPPARAM = BPPARAM)

message("Done")
