
source(".Rprofile")
suppressPackageStartupMessages({
    library("logger")
    library("readr")
    library("dplyr")
    library("DESeq2")
    library("apeglm")
})
BPPARAM = BiocParallel::MulticoreParam(snakemake@threads[[1]])

################################################################################
#
# author: Christoph Kiefer
# email: christophak@bmb.sdu.dk
#
################################################################################

log_info("Importing data...")
dds <- readRDS(snakemake@input[["dds"]])
biomart <- readRDS(snakemake@input[["biomart"]])

# set reference
colData(dds)$condition_temp <- relevel(colData(dds)$condition_temp, ref = "22")

# filtering lowly expressed genes
dds <- dds[rowSums(counts(dds)) > 10, ]

log_info("Fitting models...")
dds <- dds %>%
    DESeq(parallel = TRUE, BPPARAM = BPPARAM)

log_info("Calculating results tables...")
# results table contains s value and log2FC for the coef
res <- dds %>%
    lfcShrink(.,
        coef = "condition_temp_4_vs_22",
        lfcThreshold = snakemake@params[["lfcThreshold"]],
        type = "apeglm",
        parallel = TRUE, BPPARAM = BPPARAM)


log_info("Writing to disc...")
id <- if_else(snakemake@params$type == "gene",
        "ensembl_gene_id_version", "ensembl_transcript_id_version")

res %>%
    as_tibble(rownames = id) %>%
    left_join(biomart, by = id) %>%
    dplyr::select(starts_with("ensembl_"), mgi_symbol, description,
        gene_biotype, everything()) %>%
    write_csv(snakemake@output[[1]])


log_success("Done.")
