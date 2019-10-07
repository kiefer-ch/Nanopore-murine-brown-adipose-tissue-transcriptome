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

output_name <- args[1]
filter_by <- args[2]


# prepare sample_info
sample_info <- read_csv("sample_info/sampleInfo.csv") %>%
    filter(!is.na(get(filter_by))) %>%
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

saveRDS(cts, output_name)
