#!/usr/bin/Rscript --no-restore --no-environ --no-save

# set libpaths to packrat local library
source("packrat/init.R")
source("packrat/init.R")
library("dplyr")
library("readr")
library("purrr")

################################################################################
#
# author: Christoph Kiefer
# email: christophak@bmb.sdu.dk
#
################################################################################

# biomart
biomaRt_tx <- read_rds(snakemake@input[["biomaRt_tx"]])
biomaRt_gene <- read_rds(snakemake@input[["biomaRt_gene"]])
biotype_groups <- read_csv(snakemake@input[["biotype_groups"]])

# gene level
df_gene <- c(snakemake@input[["gene_counts"]], snakemake@input[["gene_tpm"]]) %>%
    set_names(basename(.)) %>%
    map(read_csv) %>%
    map(select, -mgi_symbol, -description, -gene_biotype) %>%
    map(tidyr::gather, key = "sample_id", value = "counts", -ensembl_gene_id_version) %>%
    bind_rows(.id = "dataset") %>%
    tidyr::separate(dataset, "dataset", extra = "drop")

df_gene_average <- df_gene %>%
    group_by(dataset, ensembl_gene_id_version) %>%
    summarise(average_counts = mean(counts)) %>%
    left_join(biomaRt_gene, by = "ensembl_gene_id_version") %>%
    left_join(biotype_groups, by = c("gene_biotype" = "transcript_biotype"))

df_gene_average %>%
    ungroup() %>%
    tidyr::spread(dataset, average_counts) %>%
    write_csv(snakemake@output[["gene"]])

# tx level
df_tx <- c(snakemake@input[["tx_counts"]], snakemake@input[["tx_tpm"]]) %>%
    set_names(basename(.)) %>%
    map(read_csv) %>%
    map(select, -mgi_symbol, -description, -gene_biotype, -transcript_biotype, -transcript_length) %>%
    map(tidyr::gather, key = "sample_id", value = "counts",
        -ensembl_gene_id_version, -ensembl_transcript_id_version) %>%
    bind_rows(.id = "dataset") %>%
    tidyr::separate(dataset, "dataset", extra = "drop")

df_tx_average <- df_tx %>%
    group_by(dataset, ensembl_transcript_id_version) %>%
    summarise(average_counts = mean(counts)) %>%
    left_join(biomaRt_tx, by = "ensembl_transcript_id_version") %>%
    left_join(biotype_groups, by = "transcript_biotype")

df_tx_average %>%
    ungroup() %>%
    tidyr::spread(dataset, average_counts) %>%
    write_csv(snakemake@output[["tx"]])
