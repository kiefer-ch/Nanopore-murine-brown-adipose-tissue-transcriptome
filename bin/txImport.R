#!/usr/bin/Rscript --no-restore --no-environ --no-save

source("packrat/init.R")
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("tximport"))
suppressPackageStartupMessages(library("DESeq2"))

################################################################################
#
# author: Christoph Kiefer
# email: christophak@bmb.sdu.dk
#
################################################################################

# generate tx2gene table
tx2g <- loadDB(snakemake@input[["txdb"]])  %>%
    AnnotationDbi::select(.,
        keys =  AnnotationDbi::keys(., keytype = "GENEID"),
        keytype = "GENEID",
        columns = "TXNAME") %>%
    dplyr::select(TXNAME, GENEID)

# prepare sample_info
sample_info <- read_csv(snakemake@input[["sample_info"]]) %>%
    filter(!is.na(illumina)) %>%
    mutate_at(vars(matches("condition")), as.factor) %>%
    mutate(path = file.path("salmon", illumina, "quant.sf")) # this line should be changed!

# create and output dds
tximport(files = sample_info$path,
        type = "salmon",
        tx2gene = tx2g,
        txOut = txout) %>%
    DESeqDataSetFromTximport(.,
        colData = tibble::column_to_rownames(as.data.frame(sample_info), "illumina"),
        design = ~ condition_temp) %>%
    saveRDS(snakemake@output[[1]])
