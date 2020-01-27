#!/usr/bin/Rscript --no-restore --no-environ --no-save

source("packrat/init.R")
library("dplyr")
library("readr")
library("DRIMSeq")
library("stageR")
BPPARAM = BiocParallel::MulticoreParam(snakemake@threads[[1]])

################################################################################
#
# author: Christoph Kiefer
# email: christophak@bmb.sdu.dk
#
################################################################################

# import data
dmds <- read_rds(snakemake@input[["dmds"]])
biomart_gene <- read_rds(snakemake@input[["biomaRt_gene"]])

# total number of samples
n <- 6
# number of samples in the smallest group
n.small <- 3

# filter
dmds <- dmFilter(dmds,
    min_samps_feature_expr = n.small, min_feature_expr = 10,
    min_samps_feature_prop = n.small, min_feature_prop = 0.1,
    min_samps_gene_expr = n, min_gene_expr = 10)

# fitAndTestModels
design_full <- model.matrix(~condition_temp,
    data = DRIMSeq::samples(dmds) %>%
        mutate(condition_temp = relevel(condition_temp, "22")))

dmds <- dmPrecision(dmds, design = design_full,
    BPPARAM = BPPARAM)
dmds <- dmFit(dmds, design = design_full)
dmds <- dmTest(dmds, coef = "condition_temp4")

# export
saveRDS(dmds, snakemake@output[["dmds"]])

# build results tables
res_gene <- DRIMSeq::results(dmds)
res_tx <- DRIMSeq::results(dmds, level = "feature")

# exchangeNAfor1
no.na <- function(x) if_else(is.na(x), 1, x)

res_gene$pvalue <- no.na(res_gene$pvalue)
res_tx$pvalue <- no.na(res_tx$pvalue)

# DrimSeq does not allow ':' in the names :-(
pScreen <- res_gene$pvalue
names(pScreen) <- sub(':', '_', res_gene$gene_id)

pConfirmation <- matrix(res_tx$pvalue, ncol = 1)
rownames(pConfirmation) <- res_tx$feature_id

tx2gene <- res_tx[ , c("feature_id", "gene_id")] %>%
    mutate(gene_id = sub(':', '_', gene_id))

stageRObj <- stageRTx(pScreen = pScreen, pConfirmation = pConfirmation,
        pScreenAdjusted = FALSE, tx2gene = tx2gene)
stageRObj <- stageWiseAdjustment(stageRObj, method = "dtu", alpha = 0.05)

drim.padj <- getAdjustedPValues(stageRObj, order = TRUE,
    onlySignificantGenes = FALSE)

res_stageR <- drim.padj %>%
    as_tibble() %>%
    dplyr::rename(ensembl_gene_id_version = "geneID",
        ensembl_transcript_id_version = "txID") %>%
    left_join(biomart_gene, by = "ensembl_gene_id_version")

# Undo the name change
res_stageR <- res_stageR %>%
    mutate(ensembl_gene_id_version = sub('_', ':', ensembl_gene_id_version))

# export
res_stageR %>%
    write_csv(snakemake@output[["res"]])
