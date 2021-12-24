
source(".Rprofile")
suppressPackageStartupMessages({
    library("logger")
    library("dplyr")
    library("readr")
    library("DRIMSeq")
    library("stageR")
})

BPPARAM = BiocParallel::MulticoreParam(snakemake@threads[[1]])

################################################################################
#
# author: Christoph Kiefer
# email: christophak@bmb.sdu.dk
#
################################################################################

log_info("Importing data...")
dmds <- read_rds(snakemake@input[["dmds"]])

sample_info <- read_csv(snakemake@input[["sample_info"]],
        show_col_types = FALSE) %>%
    mutate_at(vars(matches("condition")), as.factor) %>%
    filter(!is.na(ont)) %>%
    dplyr::select(condition_temp, sample_id = "cdna")

log_info("Filtering lowly expressed transcripts...")
genes <- counts(dmds)$gene_id %>%  unique() %>%  length()
transcripts <- counts(dmds)$feature_id %>%  unique() %>%  length()
log_info(sprintf("Before filtering: %s genes with %s transcripts", genes, transcripts))

# total number of samples
n <- 6
# number of samples in the smallest group
n.small <- 3

# filter
dmds <- dmFilter(dmds,
    min_samps_feature_expr = n.small,
    min_feature_expr = snakemake@params$min_feature_expr,
    min_samps_feature_prop = n.small,
    min_feature_prop = snakemake@params$min_feature_prop,
    min_samps_gene_expr = n,
    min_gene_expr = snakemake@params$min_gene_expr)

genes <- counts(dmds)$gene_id %>%  unique() %>%  length()
transcripts <- counts(dmds)$feature_id %>%  unique() %>%  length()
log_info(sprintf("After filtering: %s genes with %s transcripts", genes, transcripts))


log_info("Fitting and testing models...")
design_full <- model.matrix(~condition_temp,
    data = DRIMSeq::samples(dmds) %>%
        mutate(condition_temp = relevel(condition_temp, "22")))

dmds <- dmPrecision(dmds, design = design_full,
    BPPARAM = BPPARAM)
dmds <- dmFit(dmds, design = design_full)
dmds <- dmTest(dmds, coef = "condition_temp4")

# export
log_info("Writing dmds to disc...")
saveRDS(dmds, snakemake@output[["dmds"]])

log_success("Done")
