#!/usr/bin/Rscript --no-restore --no-environ --no-save

rm(list = ls())

suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("tximport"))
suppressPackageStartupMessages(library("DESeq2"))

################################################################################
#
# author: Christoph Kiefer
# email: christophak@bmb.sdu.dk
#
# usage: tximport.R sample_info annotation.gtf output.rds txout filter_by
# txout must be either TRUE or FALSE
#
################################################################################
# read command line args
args = commandArgs(trailingOnly = TRUE)

if (length(args) != 5) {
    stop("This script needs exactly 5 arguments to be called.\n")
}

sample_df <- args[1]
annotation_file <- args[2]
output_name <- args[3]

if (args[4] == "--genelevel") {
    txout <- FALSE
} else if (args[4] == "--txlevel") {
    txout <- TRUE
} else {
    stop("4th parameter must be one of '--txlevel' or '--genelevel'.")
}

filter_by <- args[5]

# generate tx2gene table
tx2g <- GenomicFeatures::makeTxDbFromGFF(annotation_file,
        format = "gtf", circ_seqs = character()) %>%
    AnnotationDbi::select(.,
        keys =  AnnotationDbi::keys(., keytype = "GENEID"),
        keytype = "GENEID",
        columns = "TXNAME") %>%
    dplyr::select(TXNAME, GENEID)

# prepare sample_info
sample_info <- read_csv(sample_df) %>%
    filter(!is.na(get(filter_by))) %>%
    mutate_at(vars(matches("condition")), as.factor) %>%
    mutate(path = file.path("salmon", illumina, "quant.sf"))

# create and output dds
tximport(files = sample_info$path,
        type = "salmon",
        tx2gene = tx2g,
        txOut = txout) %>%
    DESeqDataSetFromTximport(.,
        colData = tibble::column_to_rownames(as.data.frame(sample_info), "illumina"),
        design = ~ condition_temp) %>%
    saveRDS(output_name)
