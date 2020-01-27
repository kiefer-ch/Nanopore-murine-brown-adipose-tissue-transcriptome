#!/usr/bin/Rscript --no-restore --no-environ --no-save

source("packrat/init.R")
library("readr")
library("dplyr")

################################################################################
#
# author: Christoph Kiefer
# email: christophak@bmb.sdu.dk
#
################################################################################

message("Preparing tx2g...")
tx2g <- AnnotationDbi::loadDb(snakemake@input[["txdb"]]) %>%
    AnnotationDbi::select(.,
        keys =  AnnotationDbi::keys(., keytype = "GENEID"),
        keytype = "GENEID",
        columns = "TXNAME") %>%
    dplyr::select(TXNAME, GENEID) %>%
    dplyr::rename(ensembl_gene_id_version = "GENEID",
        ensembl_transcript_id_version = "TXNAME")

message("Collecting sample info")
sample_info <- suppressMessages(read_csv(snakemake@input[["sampleInfo"]])) %>%
    mutate_at(vars(matches("condition")), as.factor) %>%
    dplyr::select(sample_id, condition_temp)

message("Calculating which transcripts should be kept...")
# keep only those transcripts, that in at least one condition have more then
# threshold % of the TPM
df <- read_csv(snakemake@input[["counts"]]) %>%
    dplyr::select(-mgi_symbol, -ensembl_gene_id_version, -description,
        -gene_biotype, -transcript_biotype, -transcript_length) %>%
    left_join(tx2g, by = "ensembl_transcript_id_version") %>%
    tidyr::gather(key = sample_id, value = "counts", -ensembl_gene_id_version,
        -ensembl_transcript_id_version) %>%
    mutate(counts = (2^counts) - 1) %>% # ntd is log2(x+1)
    filter(counts != 0) %>%
    left_join(sample_info, by = "sample_id") %>%
    group_by(condition_temp, ensembl_gene_id_version, ensembl_transcript_id_version) %>%
    summarise(counts = sum(counts)) %>%
    mutate(freq = counts / sum(counts)) %>%
    filter(freq > snakemake@params[["threshold"]] / 100)

keep_tx <- df %>%
    pull(ensembl_transcript_id_version) %>%
    unique()

message("Reading gtf file...")
read_gtf <- function(file) {
    read_delim(file = file,
        delim = "\t",
        comment = "#",
        na = c('.'),
        col_names = c("sequence", "source", "feature", "start", "end", "score",
            "strand", "phase", "attributes"),
        col_types = cols(
            sequence = col_character(),
            source = col_character(),
            feature = col_character(),
            start = col_integer(),
            end = col_integer(),
            score = col_character(),
            strand = col_character(),
            phase = col_character(),
            attributes = col_character()),
        progress = FALSE)
}

gtf <- read_gtf(snakemake@input[["annotation"]])

message("Subsetting gtf...")
gtf <- gtf %>%
    tidyr::separate(attributes, c("gene_id", "transcript_id"),
        sep = ';', remove = FALSE, extra = "drop") %>%
    mutate(transcript_id = if_else(feature == "gene", NA_character_, substr(transcript_id, 17, nchar(transcript_id) - 1))) %>%
    filter(feature == "gene" | transcript_id %in% keep_tx)

keep_gene <- gtf %>%
    filter(feature == "transcript") %>%
    pull(gene_id)

gtf <- gtf %>%
    filter(gene_id %in% keep_gene) %>%
    dplyr::select(-gene_id, -transcript_id)

message("Writing to disk...")
withr::with_options(c(scipen = 10),
    write.table(gtf, snakemake@output[[1]],
        col.names = FALSE, row.names = FALSE,
        na = '.', quote = FALSE, sep = '\t'))

message("Done")
