#!/usr/bin/Rscript --no-restore --no-environ --no-save

# set libpaths to packrat local library
source("packrat/init.R")
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("DEXSeq"))
suppressPackageStartupMessages(library("readr"))
    select <- dplyr::select

################################################################################
#
# author: Christoph Kiefer
# email: christophak@bmb.sdu.dk
#
################################################################################

dxd <- read_rds(snakemake@input[["dxd"]])
dxr <- DEXSeqResults(dxd)

biomart <- read_rds(snakemake@input[["biomaRt_gene"]])

dxr %>%
    as_tibble() %>%
    group_by(groupID) %>%
    summarise(pvalue = min(pvalue),
        padj = min(padj)) %>%
    arrange(pvalue) %>%
    dplyr::rename(ensembl_gene_id_version = "groupID") %>%
    left_join(biomart, by = "ensembl_gene_id_version") %>%
    select(ensembl_gene_id_version, mgi_symbol, everything()) %>%
    write_csv(snakemake@output[[1]])
