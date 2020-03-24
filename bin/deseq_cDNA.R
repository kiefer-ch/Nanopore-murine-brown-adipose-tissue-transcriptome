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
sample_info <- read_csv(snakemake@input[["sample_info"]]) %>%
    mutate_at(vars(matches("condition")), as.factor) %>%
    filter(!is.na(cdna))

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
    tidyr::separate(file, sep = '_', c("date", "pool", "waste", "barcode"),
        extra = "drop") %>%
    select(-date, -waste) %>%
    group_by(barcode, ensembl_gene_id_version, ensembl_transcript_id_version) %>%
    summarise(count = sum(count)) %>%
    tidyr::spread(key = barcode, value = count) %>%
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

# set colnames to sample_id
lookup <- sample_info$sample_id
names(lookup) <- sample_info$cdna

colnames(counts) <- unname(lookup[colnames(counts)])

# create dds
dds <- DESeqDataSetFromMatrix(countData = counts,
    colData = tibble::column_to_rownames(as.data.frame(sample_info), "sample_id"),
    design = as.formula(snakemake@params[["design"]]))
saveRDS(dds, snakemake@output[["dds"]])
