#!/usr/bin/Rscript --no-restore --no-environ --no-save

# set libpaths to packrat local library
save.image()
source(".Rprofile")
library("dplyr")
library("DEXSeq")
library("readr")
library("stageR")
    select <- dplyr::select

################################################################################
#
# author: Christoph Kiefer
# email: christophak@bmb.sdu.dk
#
################################################################################

# import data
dxd <- read_rds(snakemake@input[["dxd"]])
biomart <- read_rds(snakemake@input[["biomaRt_gene"]])

# DEXSeq results
dxr <- DEXSeqResults(dxd)
qval <- perGeneQValue(dxr)

# stageR
pScreen <- qval
pConfirmation <- matrix(dxr$pvalue, ncol = 1)
dimnames(pConfirmation) <- list(rownames(dxr), "exon")
# stageR does not allow ':' in the names :-(
rownames(pConfirmation) <- sub(':', '_', rownames(pConfirmation))

pScreen <- tidyr::replace_na(pScreen, 1)
pConfirmation <- tidyr::replace_na(pConfirmation, 1)

tx2gene <- as.data.frame(dxr[,c("featureID", "groupID")]) %>%
    select(-featureID) %>%
    tibble::rownames_to_column("featureID") %>%
    mutate(featureID = sub(':', '_', featureID))

stageRObj <- stageRTx(pScreen = pScreen, pConfirmation = pConfirmation,
                      pScreenAdjusted = TRUE, tx2gene = tx2gene)
stageRObj <- stageWiseAdjustment(stageRObj, method = "dtu", alpha = 0.05)

dex.padj <- getAdjustedPValues(stageRObj, order = TRUE,
                                onlySignificantGenes = FALSE)

# data output
res_stageR <- dex.padj %>%
    as_tibble() %>%
    dplyr::rename(ensembl_gene_id_version = "geneID",
                  exon_number = "txID") %>%
    left_join(biomart, by = "ensembl_gene_id_version")

res_stageR <- res_stageR %>%
    mutate(ensembl_gene_id_version = sub('_', ':', ensembl_gene_id_version))

res_stageR %>%
    write_csv(snakemake@output[[1]])
