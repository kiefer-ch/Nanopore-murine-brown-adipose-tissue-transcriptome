#!/usr/bin/Rscript --no-restore --no-environ --no-save

# set libpaths to packrat local library
source("../packrat/init.R")

suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("dplyr"))

################################################################################
#
# author: Christoph Kiefer
# email: christophak@bmb.sdu.dk
#
################################################################################

df <- read_rds(snakemake@input[["scaledTPM"]]) %>%
    as_tibble(rownames = "ensembl_transcript_id_version")

tx2g <- GenomicFeatures::makeTxDbFromGFF(snakemake@input[["annotation"]],
        format = "gtf", circ_seqs = character()) %>%
    AnnotationDbi::select(.,
        keys =  AnnotationDbi::keys(., keytype = "GENEID"),
        keytype = "GENEID",
        columns = "TXNAME") %>%
    dplyr::select(TXNAME, GENEID) %>%
    dplyr::rename(ensembl_gene_id_version = "GENEID",
        ensembl_transcript_id_version = "TXNAME")

df <- df %>%
    left_join(tx2g, by = "ensembl_transcript_id_version") %>%
    tidyr::gather(key = sample, value = "scaledTPM", -ensembl_gene_id_version,
        -ensembl_transcript_id_version) %>%
    filter(scaledTPM != 0) %>%
    group_by(sample, ensembl_gene_id_version) %>%
    mutate(freq = scaledTPM / sum(scaledTPM)) %>%
    filter(freq > .15)

keep_tx <- df %>%
    pull(ensembl_transcript_id_version) %>%
    unique()

read_gtf <- function(file) {
    vroom::vroom(file = file,
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
        num_threads = snakemake@threads[[1]],
        progress = FALSE)
}

gtf <- read_gtf(snakemake@input[["annotation"]])

gtf <- gtf[1:100,] %>%
    tidyr::separate(attributes, c("gene_id", "transcript_id"),
        sep = ';', remove = FALSE, extra = "drop") %>%
    mutate(transcript_id = if_else(feature == "gene", NA_character_, substr(transcript_id, 17, nchar(transcript_id) - 1))) %>%
    filter(feature == "gene" | transcript_id %in% keep_tx)

keep_gene <- gtf %>%
    filter(feature == "transcript") %>%
    pull(gene_id)

gtf %>%
    filter(gene_id %in% keep_gene) %>%
    dplyr::select(-gene_id, -transcript_id) %>%
    write_tsv(snakemake@output[[1]],
        na = '.', col_names = FALSE)
