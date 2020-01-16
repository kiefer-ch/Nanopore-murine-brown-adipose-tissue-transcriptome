#!/usr/bin/Rscript --no-restore --no-environ --no-save

source("packrat/init.R")
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("DESeq2"))
    select <- dplyr::select
    rename <- dplyr::rename
BPPARAM = BiocParallel::MulticoreParam(snakemake@threads[[1]])

################################################################################
#
# author: Christoph Kiefer
# email: christophak@bmb.sdu.dk
#
################################################################################

# prepare sample_info
sample_info <- read_csv(snakemake@input[["sample_info"]]) %>%
    mutate_at(vars(matches("condition")), as.factor) %>%
    filter(!is.na(ont))

# import data
counts <- read_tsv(snakemake@input[["counts"]]) %>%
    tidyr::separate(transcript, c("ensembl_transcript_id_version", "ensembl_gene_id_version"),
        sep = '\\|', extra = "drop") %>%
    tidyr::drop_na()

# tx level
if (as.logical(snakemake@params[["txOut"]])) {
    counts <- counts %>%
        select(-ensembl_gene_id_version) %>%
        as.data.frame() %>%
        tibble::column_to_rownames("ensembl_transcript_id_version") %>%
        data.matrix()
} else if (!as.logical(snakemake@params[["txOut"]])) {
# gene level
    counts <- counts %>%
        group_by(ensembl_gene_id_version) %>%
        summarise_if(is.numeric, sum) %>%
        tibble::column_to_rownames("ensembl_gene_id_version") %>%
        data.matrix()
}

lookup <- sample_info$sample_id
names(lookup) <- sample_info$ont

colnames(counts) <- unname(lookup[colnames(counts)])

# create dds
dds <- DESeqDataSetFromMatrix(countData = counts,
    colData = tibble::column_to_rownames(as.data.frame(sample_info), "sample_id"),
    design = as.formula(snakemake@params[["design"]]))
saveRDS(dds, snakemake@output[["dds"]])
