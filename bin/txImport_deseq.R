#!/usr/bin/Rscript --no-restore --no-environ --no-save

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

save.image("tximport.RData")

# generate tx2gene table
tx2g <- AnnotationDbi::loadDb(snakemake@input[["txdb"]])  %>%
    AnnotationDbi::select(.,
        keys =  AnnotationDbi::keys(., keytype = "GENEID"),
        keytype = "GENEID",
        columns = "TXNAME") %>%
    dplyr::select(TXNAME, GENEID)

# prepare sample_info
sample_info <- read_csv(snakemake@input[["sample_info"]]) %>%
    mutate_at(vars(matches("condition")), as.factor) %>%
    mutate(path = file.path("salmon", sample, "quant.sf")) # this line should be changed!

# create and output dds
tximport(files = sample_info$path,
        type = "salmon",
        tx2gene = tx2g,
        txOut = as.logical(snakemake@params[["txOut"]])) %>%
    DESeqDataSetFromTximport(.,
        colData = tibble::column_to_rownames(as.data.frame(sample_info), "sample"),
        design = as.formula(snakemake@params[["design"]])) %>%
    saveRDS(snakemake@output[[1]])
