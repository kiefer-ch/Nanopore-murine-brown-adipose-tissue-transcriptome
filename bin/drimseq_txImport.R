#!/usr/bin/Rscript --no-restore --no-environ --no-save

source("packrat/init.R")
library("readr")
library("dplyr")
library("tximport")

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

# generate tx2gene table
tx2g <- AnnotationDbi::loadDb(snakemake@input[["txdb"]])  %>%
    AnnotationDbi::select(.,
        keys =  AnnotationDbi::keys(., keytype = "GENEID"),
        keytype = "GENEID",
        columns = "TXNAME") %>%
    dplyr::select(TXNAME, GENEID)

# create and output dds
txi <- tximport(files = sample_info$path,
    type = "salmon",
    tx2gene = tx2g,
    txOut = TRUE,
    countsFromAbundance = "dtuScaledTPM")

cts <- txi$counts
cts <- cts[rowSums(cts) > 0, ]
colnames(cts) <- sample_info$sample_id

saveRDS(cts, snakemake@output[[1]])
