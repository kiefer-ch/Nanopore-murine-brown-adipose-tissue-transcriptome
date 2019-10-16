#!/usr/bin/Rscript --no-restore --no-environ --no-save

# set libpaths to packrat local library
source("packrat/init.R")

suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("dplyr"))

################################################################################
#
# author: Christoph Kiefer
# email: christophak@bmb.sdu.dk
#
################################################################################

message("Preparing tx2g...")
tx2g <- read_rds(path = snakemake@input[["txdb"]]) %>%
    dplyr::select(TXNAME, GENEID) %>%
    dplyr::rename(ensembl_gene_id_version = "GENEID",
        ensembl_transcript_id_version = "TXNAME")

message("Calculating which transcripts should be kept...")
df <- read_rds(snakemake@input[["scaledTPM"]]) %>%
    as_tibble(rownames = "ensembl_transcript_id_version") %>%
    left_join(tx2g, by = "ensembl_transcript_id_version") %>%
    tidyr::gather(key = sample, value = "scaledTPM", -ensembl_gene_id_version,
        -ensembl_transcript_id_version) %>%
    filter(scaledTPM != 0) %>%
    group_by(sample, ensembl_gene_id_version) %>%
    mutate(freq = scaledTPM / sum(scaledTPM)) %>%
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
