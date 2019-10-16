#!/usr/bin/Rscript --no-restore --no-environ --no-save

# set libpaths to packrat local library
source("../packrat/init.R")

suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("DEXSeq"))

################################################################################
#
# author: Christoph Kiefer
# email: christophak@bmb.sdu.dk
#
################################################################################
sample_info <- read_csv(snakemake@input[["sample_info"]]) %>%
    filter(!is.na(illumina)) %>%
    mutate_at(vars(matches("condition")), as.factor) %>%
    mutate(path = file.path("dexseq", paste0(illumina, ".txt")))

DEXSeqDataSetFromHTSeq(
    countfiles = sample_info$path,
    sampleData = tibble::column_to_rownames(as.data.frame(sample_info), "illumina"),
    design = ~ condition_temp + exon + condition_temp:exon,
    flattenedfile = snakemake@input[["annotation"]]) %>%
    saveRDS(snakemake@output[[1]])
