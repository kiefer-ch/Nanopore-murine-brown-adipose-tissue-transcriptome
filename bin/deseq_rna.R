#!/usr/bin/Rscript --no-restore --no-environ --no-save

source(".Rprofile")
library("dplyr")
library("readr")
library("DESeq2")
library("purrr")
    select <- dplyr::select
    rename <- dplyr::rename

################################################################################
#
# author: Christoph Kiefer
# email: christophak@bmb.sdu.dk
#
################################################################################

# prepare sample_info
sample_info <- tibble(sample_id = c("rt", "cool"), condition_temp = c("22", "4"))

# import data
counts <- snakemake@input[["counts"]] %>%
    set_names(.) %>%
    map(read_tsv, col_names = FALSE) %>%
    map(select_at, vars(1:2, 11)) %>%
    map(tidyr::drop_na) %>%
    map(rename, ensembl_gene_id_version = "X2",
        ensembl_transcript_id_version = "X1",
        count = "X11") %>%
    bind_rows(.id = "file") %>%
    mutate(file = basename(file)) %>%
    tidyr::separate(file, sep = '_', c("waste", "temperature"),
        extra = "drop") %>%
    select(-waste) %>%
    group_by(temperature, ensembl_gene_id_version, ensembl_transcript_id_version) %>%
    summarise(count = sum(count)) %>%
    tidyr::spread(key = temperature, value = count) %>%
    ungroup()

# tx level
if (as.logical(snakemake@params[["txOut"]])) {
    counts <- counts %>%
        select(-ensembl_gene_id_version) %>%
        as.data.frame() %>%
        tibble::column_to_rownames("ensembl_transcript_id_version") %>%
        data.matrix()
} else if (!as.logical(snakemake@params[["txOut"]])) {
# gene level
    counts <- counts %>%
        group_by(ensembl_gene_id_version) %>%
        summarise_if(is.numeric, sum) %>%
        tibble::column_to_rownames("ensembl_gene_id_version") %>%
        data.matrix()
}

# create dds
dds <- DESeqDataSetFromMatrix(countData = counts,
    colData = sample_info,
    design = as.formula(snakemake@params[["design"]]))
saveRDS(dds, snakemake@output[["dds"]])
