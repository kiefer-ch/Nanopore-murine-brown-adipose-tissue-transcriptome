#!/usr/bin/Rscript --no-restore --no-environ --no-save

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
sample_info <- read_csv(snakemake@input[["sample_info"]]) %>%
    filter(!is.na(ont)) %>%
    mutate_at(vars(matches("condition")), as.factor) %>%
    mutate(path = file.path("salmon", illumina, "quant.sf"))

# create and output dds
txi <- tximport(files = sample_info$path,
    type = "salmon",
    txOut = TRUE,
    countsFromAbundance = "scaledTPM")

cts <- txi$counts
cts <- cts[rowSums(cts) > 0, ]
colnames(cts) <- sample_info$sample_id

saveRDS(cts, snakemake@output[[1]])
