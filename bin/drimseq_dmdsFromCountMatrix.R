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
cts <- read_csv(snakemake@input[["counts"]]) %>%
    dplyr::select(-ensembl_gene_id_version, -mgi_symbol, -description, -gene_biotype, -transcript_biotype,
        -transcript_length) %>%
    dplyr::rename(TXNAME = "ensembl_transcript_id_version")

txdf <- txdf[match(cts$TXNAME, txdf$TXNAME),]

counts <- cts %>%
    left_join(txdf, by = "TXNAME") %>%
    dplyr::rename(gene_id = "GENEID",
        feature_id = "TXNAME") %>%
    as.data.frame()

message("Import sample info...")
sample_info <- read_csv(snakemake@input[["sample_info"]]) %>%
    filter(!is.na(ont)) %>%
    mutate_at(vars(matches("condition")), as.factor) %>%
    as.data.frame()

message("Create dmds")
dmds <- dmDSdata(counts = counts, samples = sample_info)
saveRDS(dmds, snakemake@output[[1]])
