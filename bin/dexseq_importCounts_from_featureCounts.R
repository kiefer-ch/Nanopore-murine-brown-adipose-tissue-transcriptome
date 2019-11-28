#!/usr/bin/Rscript --no-restore --no-environ --no-save

# set libpaths to packrat local library
source("packrat/init.R")
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("DEXSeq"))
source("bin/load_SubreadOutput.R")

################################################################################
#
# author: Christoph Kiefer
# email: christophak@bmb.sdu.dk
#
################################################################################
message("Collecting sample info...")
sample_info <- read_csv(snakemake@input[["sample_info"]]) %>%
    filter(!is.na(ont)) %>%
    mutate_at(vars(matches("condition")), as.factor) %>%
    as.data.frame() %>%
    tibble::column_to_rownames("illumina")

message("Creating DEXSeqDataSet...")
DEXSeqDataSetFromFeatureCounts(snakemake@input[["counts"]],
    flattenedfile = snakemake@input[["annotation"]],
    sampleData = sample_info,
    design = ~ sample + exon + condition_temp:exon) %>%
    saveRDS(snakemake@output[[1]])

message("Done")
