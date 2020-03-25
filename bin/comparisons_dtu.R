#!/usr/bin/Rscript --no-restore --no-environ --no-save

# set libpaths to packrat local library
save.image()
source(".Rprofile")
library("dplyr")
library("readr")
library("purrr")

################################################################################
#
# author: Christoph Kiefer
# email: christophak@bmb.sdu.dk
#
################################################################################
# dexseq
dexseq <- snakemake@input[["dexseq"]] %>%
    set_names(., tools::file_path_sans_ext(basename(.), compression = TRUE)) %>%
    map(read_csv) %>%
    map(select, -mgi_symbol, -transcript, -description) %>%
    bind_rows(.id = "dataset") %>%
    group_by(dataset, ensembl_gene_id_version) %>%
    summarise(gene_biotype = unique(gene_biotype),
        padj = min(gene)) %>%
    tidyr::separate(dataset, c("dataset", "library"), extra = "drop") %>%
    tidyr::unite("method", dataset, library) %>%
    tidyr::spread(method, padj)

# drimseq
drimseq <- snakemake@input[["drimseq"]] %>%
    set_names(., tools::file_path_sans_ext(basename(.), compression = TRUE)) %>%
    map(read_csv) %>%
    map(select, -mgi_symbol, -transcript, -description) %>%
    bind_rows(.id = "dataset") %>%
    group_by(dataset, ensembl_gene_id_version) %>%
    summarise(gene_biotype = unique(gene_biotype),
        padj = min(gene)) %>%
    tidyr::separate(dataset, c("dataset", "flair"), extra = "drop") %>%
    mutate(library = "drimseq") %>%
    mutate(dataset = if_else(flair == "flair", paste0(dataset, "_flair"), dataset)) %>%
    select(-flair) %>%
    tidyr::unite("method", dataset, library) %>%
    tidyr::spread(method, padj) %>%
    filter(substr(ensembl_gene_id_version, 1, 3) == "ENS")

df <- list(drimseq, dexseq) %>%
    reduce(full_join, by = c("ensembl_gene_id_version", "gene_biotype")) %>%
    mutate_if(is.numeric, function(x) tidyr::replace_na(x, 1))

df %>%
    saveRDS(snakemake@output[["all"]])

biomart <- read_rds(snakemake@input[["biomaRt_gene"]])

df_sub <- df %>%
    select(-gene_biotype) %>%
    left_join(biomart, by = "ensembl_gene_id_version") %>%
    select(ensembl_gene_id_version, mgi_symbol, description, gene_biotype, everything()) %>%
    filter_if(is.numeric, any_vars(. < .05))

counts <- read_csv(snakemake@input[["counts"]]) %>%
    dplyr::select(-mgi_symbol, -description, -gene_biotype, -grouped_biotype) %>%
    tidyr::gather("dataset", "counts", -ensembl_gene_id_version) %>%
    mutate(counts = 2^counts - 1) %>%
    mutate(dataset = if_else(dataset == "illumina",
        paste0("avgTpm_", dataset),
        paste0("avgCounts_", dataset))) %>%
    tidyr::spread(dataset, counts)

df_sub %>%
    rename_if(is.numeric, paste0, "_padj") %>%
    left_join(counts, by = "ensembl_gene_id_version") %>%
    dplyr::select(ensembl_gene_id_version, mgi_symbol, description,
        gene_biotype, avgCounts_cdna, avgCounts_teloprime,
        avgTpm_illumina, everything()) %>%
    mutate(sort = mean(c(avgCounts_cdna, avgCounts_teloprime,
                       avgTpm_illumina))) %>%
    arrange(desc(sort)) %>%
    select(-sort) %>%
    write_csv(snakemake@output[["signif"]])
