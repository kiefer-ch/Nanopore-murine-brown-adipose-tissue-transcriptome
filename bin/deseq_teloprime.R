#!/usr/bin/Rscript --no-restore --no-environ --no-save

# set libpaths to packrat local library
message("Load packages...")
source("packrat/init.R")
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("DESeq2"))
suppressPackageStartupMessages(library("Mus.musculus"))
    select <- dplyr::select
    rename <- dplyr::rename

################################################################################
#
# author: Christoph Kiefer
# email: christophak@bmb.sdu.dk
#
################################################################################

save.image("deseq_teloprime.RData")

# import data
message("Import data...")
raw <- read_tsv(snakemake@input[["ont_gene_raw"]])
sample_info <- read_csv(snakemake@input[["sample_info"]]) %>%
    filter(!is.na(ont)) %>%
    mutate_at(vars(matches("condition")), as.factor)

# tidy rownames
raw <- raw%>%
    tidyr::separate(transcript, c("transcript_id_ens", "gene_id_ens"), sep = '\\|') %>%
    tidyr::drop_na()

# tidy colnames
lookup <- sample_info$illumina
names(lookup) <- sample_info$ont

colnames(raw) <- lookup[colnames(raw)]
colnames(raw)[c(1, 2)] <- c("transcript_id_ens", "gene_id_ens")

# coerce to numeric matrix
raw_tx <- raw %>%
    select(-gene_id_ens) %>%
    as.data.frame() %>%
    magrittr::set_rownames(.$transcript_id_ens) %>%
    `[`( ,-1) %>%
    data.matrix()

# create dds
dds <- DESeqDataSetFromMatrix(countData = raw_tx,
    colData = sample_info,
    design = ~ condition_temp)

################################################################################
# txlevel
################################################################################

message("Export txlevel count matrices...")

# prepare countmatrices
rld <- rlog(dds, blind = FALSE)
ntd <- normTransform(dds)

# annotate count matrices
annotate_cm <- function(tibble, keytype = "ENSEMBL") {
    ks <- tibble$gene_id_ens %>%
        tools::file_path_sans_ext()

    sym <- mapIds(org.Mm.eg.db, keys = ks, keytype = keytype,
        column = "SYMBOL", multiVals = "first") %>%
        tibble::enframe(name = "gene_id", value = "mgi_symbol")
    names <- mapIds(org.Mm.eg.db, keys = ks, keytype = keytype,
        column = "GENENAME", multiVals = "first") %>%
        tibble::enframe(name = "gene_id", value = "gene_name")

    tibble %>%
        mutate(gene_id = tools::file_path_sans_ext(gene_id_ens)) %>%
        left_join(sym, by = "gene_id") %>%
        left_join(names, by = "gene_id") %>%
        select(-gene_id) %>%
        select(gene_id_ens, mgi_symbol, gene_name, everything())
}

# export count matrices
# with variance shrinking
assay(rld) %>%
    as_tibble(rownames = "gene_id_ens") %>%
    annotate_cm(keytype = "ENSEMBLTRANS") %>%
    rename(transcript_id_ens = "gene_id_ens") %>%
    write_csv(snakemake@output[["tx_rld"]])

# without
assay(ntd) %>%
    as_tibble(rownames = "gene_id_ens") %>%
    annotate_cm(keytype = "ENSEMBLTRANS") %>%
    rename(transcript_id_ens = "gene_id_ens") %>%
    write_csv(snakemake@output[["tx_ntd"]])

counts(dds, normalized = FALSE) %>%
    as_tibble(rownames = "gene_id_ens") %>%
    annotate_cm(keytype = "ENSEMBLTRANS") %>%
    write_csv(snakemake@output[["tx_cts"]])

################################################################################
# genelevel
################################################################################

message("Export genelevel count matrices...")
# summarise tx to gene level
raw_gene <- raw %>%
        group_by(gene_id_ens) %>%
        summarise_if(is.numeric, sum) %>%
    as.data.frame() %>%
    magrittr::set_rownames(.$gene_id_ens) %>%
    `[`( ,-1) %>%
    data.matrix()

# create dds
dds <- DESeqDataSetFromMatrix(countData = raw_gene,
    colData = sample_info,
    design = ~ condition_temp)

# prepare countmatrices
rld <- rlog(dds, blind = FALSE)
ntd <- normTransform(dds)

# annotate count matrices
annotate_cm <- function(tibble) {
    ks <- tibble$gene_id_ens %>%
        tools::file_path_sans_ext()

    sym <- mapIds(org.Mm.eg.db, keys = ks, keytype = "ENSEMBL",
        column = "SYMBOL", multiVals = "first") %>%
        tibble::enframe(name = "gene_id", value = "mgi_symbol")
    names <- mapIds(org.Mm.eg.db, keys = ks, keytype = "ENSEMBL",
        column = "GENENAME", multiVals = "first") %>%
        tibble::enframe(name = "gene_id", value = "gene_name")

    tibble %>%
        mutate(gene_id = tools::file_path_sans_ext(gene_id_ens)) %>%
        left_join(sym, by = "gene_id") %>%
        left_join(names, by = "gene_id") %>%
        select(-gene_id) %>%
        select(gene_id_ens, mgi_symbol, gene_name, everything())
}

# export count matrices
# with variance shrinking
assay(rld) %>%
    as_tibble(rownames = "gene_id_ens") %>%
    annotate_cm() %>%
    write_csv(snakemake@output[["gene_rld"]])

# without
assay(ntd) %>%
    as_tibble(rownames = "gene_id_ens") %>%
    annotate_cm() %>%
    write_csv(snakemake@output[["gene_ntd"]])

# raw counts
counts(dds, normalized = FALSE) %>%
    as_tibble(rownames = "gene_id_ens") %>%
    annotate_cm() %>%
    write_csv(snakemake@output[["gene_cts"]])
