#!/usr/bin/Rscript --no-restore --no-environ --no-save

# set libpaths to packrat local library
source("packrat/init.R")
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("DEXSeq"))

save.image("dexseq_import_illumina.RData")

################################################################################
#
# author: Christoph Kiefer
# email: christophak@bmb.sdu.dk
#
################################################################################
message("Collecting sample info...")
sample_info <- read_csv(snakemake@input[["sample_info"]]) %>%
    filter(!is.na(ont)) %>%
    mutate_at(vars(matches("condition")), as.factor)

files <- snakemake@input[["files"]]
files <- files[grep(paste(pull(sample_info, snakemake@params[["dataset"]]),
    collapse = "|"), files)]

sample_info <- sample_info %>%
    mutate(path = files)

message("Creating DEXSeqDataSet...")
DEXSeqDataSetFromHTSeq(
    countfiles = sample_info$path,
    sampleData = tibble::column_to_rownames(as.data.frame(sample_info), "illumina"),
    design = ~ sample + exon + condition_temp:exon,
    flattenedfile = snakemake@input[["annotation"]]) %>%
    saveRDS(snakemake@output[[1]])

message("Done")

