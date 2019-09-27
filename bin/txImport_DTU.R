#!/usr/bin/Rscript --no-restore --no-environ --no-save

rm(list = ls())

suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("tximport"))

################################################################################
#
# author: Christoph Kiefer
# email: christophak@bmb.sdu.dk
#
################################################################################

# prepare sample_info
sample_info <- read_csv(sample_df) %>%
    filter(!is.na(get(filter_by))) %>%
    mutate_at(vars(matches("condition")), as.factor) %>%
    mutate(path = file.path("salmon", illumina, "quant.sf"))

# create and output dds
tximport(files = sample_info$path,
    type = "salmon",
    tx2gene = tx2g,
    txOut = TRUE,
    countsFromAbundance = "scaledTPM")

cts <- txi$counts
cts <- cts[rowSums(cts) > 0,]

write_csv("../data/scaledTPM.csv.gz")
