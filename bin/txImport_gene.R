#!/usr/bin/Rscript --no-restore --no-environ --no-save

rm(list = ls())

library("readr")
library("dplyr")
library("tximport")
library("DESeq2")

################################################################################
#
# author: Christoph Kiefer
# email: christophak@bmb.sdu.dk
#
################################################################################
# read command line args
args = commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
    stop("At least one argument must be supplied (sample_info.csv).\n", call. = FALSE)
}

sample_df <- args[1]
annotation_file <- args[2]
output_name <- args[3]

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
    filter(!is.na(illumina)) %>%
    mutate_at(vars(matches("condition")), as.factor) %>%
    mutate(path = file.path("salmon", illumina, "quant.sf"))

# create and output dds
tximport(files = sample_info$path, type = "salmon", tx2gene = tx2g) %>%
    DESeqDataSetFromTximport(.,
        colData = tibble::column_to_rownames(as.data.frame(sample_info), "illumina"),
        design = ~ condition_temp) %>%
    saveRDS(output_name)
