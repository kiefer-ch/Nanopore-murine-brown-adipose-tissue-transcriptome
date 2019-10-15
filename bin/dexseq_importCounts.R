#!/usr/bin/Rscript --no-restore --no-environ --no-save

.libPaths(c("packrat/lib/x86_64-pc-linux-gnu/3.6.1/",
    "packrat/lib-ext/x86_64-pc-linux-gnu/3.6.1/",
    "packrat/lib-R/x86_64-pc-linux-gnu/3.6.1/"))

suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("DEXSeq"))

################################################################################
#
# author: Christoph Kiefer
# email: christophak@bmb.sdu.dk
#
# usage: tximport.R sample_info annotation.gtf output.rds txout filter_by
# txout must be either TRUE or FALSE
#
################################################################################
sample_info <- read_csv(snakemake@input[["sample_info"]]) %>%
    filter(!is.na(get(snakemake@params[["filter_by"]]))) %>%
    mutate_at(vars(matches("condition")), as.factor) %>%
    mutate(path = file.path("dexseq", paste0(illumina, ".txt")))

DEXSeqDataSetFromHTSeq(
    countfiles = sample_info$path,
    sampleData = tibble::column_to_rownames(as.data.frame(sample_info), "illumina"),
    design = ~ condition_temp + exon + condition_temp:exon,
    flattenedfile = snakemake@input[["annotation"]]) %>%
    saveRDS(snakemake@output[[1]])
