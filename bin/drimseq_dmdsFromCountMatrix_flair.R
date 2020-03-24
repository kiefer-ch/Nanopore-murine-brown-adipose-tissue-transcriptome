#!/usr/bin/Rscript --no-restore --no-environ --no-save

# set libpaths to packrat local library
source(".Rprofile")
library("readr")
library("dplyr")
library("DRIMSeq")

################################################################################
#
# author: Christoph Kiefer
# email: christophak@bmb.sdu.dk
#
################################################################################

message("Importing txdb...")
txdf <- AnnotationDbi::loadDb(snakemake@input[["txdb"]])  %>%
    AnnotationDbi::select(.,
        keys =  AnnotationDbi::keys(., keytype = "GENEID"),
        keytype = "GENEID",
        columns = "TXNAME") %>%
    dplyr::select(TXNAME, GENEID) %>%
    group_by(GENEID) %>%
    mutate(ntx = row_number()) %>%
    ungroup()

message("Import counts...")
cts <- read_tsv(snakemake@input[["counts"]]) %>%
    tidyr::separate(ids, c("TXNAME", "GENEID"), sep = '_')
txdf <- txdf[match(cts$TXNAME, txdf$TXNAME),]

counts <- cts %>%
    left_join(txdf, by = c("TXNAME", "GENEID")) %>%
    dplyr::rename(gene_id = "GENEID",
        feature_id = "TXNAME") %>%
    as.data.frame()

message("Import sample info...")
sample_info <- read_csv(snakemake@input[["sample_info"]]) %>%
    filter(!is.na(ont)) %>%
    mutate_at(vars(matches("condition")), as.factor) %>%
    as.data.frame()

# This line is not good style!
colnames(counts)[3:8] <- sample_info$sample_id

message("Create dnds")
dmds <- dmDSdata(counts = counts, samples = sample_info)
saveRDS(dmds, snakemake@output[[1]])

