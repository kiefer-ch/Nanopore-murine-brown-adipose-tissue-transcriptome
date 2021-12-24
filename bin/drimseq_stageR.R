
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
biomart_gene <- read_rds(snakemake@input[["biomaRt_gene"]])

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

# build results tables
log_info("Building results table...")
res_gene <- DRIMSeq::results(dmds)
res_tx <- DRIMSeq::results(dmds, level = "feature")

# exchangeNAfor1
res_gene$pvalue <- tidyr::replace_na(res_gene$pvalue, 1)
res_tx$pvalue <- tidyr::replace_na(res_tx$pvalue, 1)

# stageR does not allow ':' in the names :-(
pScreen <- res_gene$pvalue
names(pScreen) <- sub(':', '_', res_gene$gene_id)

pConfirmation <- matrix(res_tx$pvalue, ncol = 1)
rownames(pConfirmation) <- res_tx$feature_id

tx2gene <- res_tx[ , c("feature_id", "gene_id")] %>%
    mutate(gene_id = sub(':', '_', gene_id))

stageRObj <- stageRTx(pScreen = pScreen, pConfirmation = pConfirmation,
        pScreenAdjusted = FALSE, tx2gene = tx2gene)
stageRObj <- stageWiseAdjustment(stageRObj, method = "dtu", alpha = 0.05)

drim.padj <- getAdjustedPValues(stageRObj, order = TRUE,
    onlySignificantGenes = FALSE)

res_stageR <- drim.padj %>%
    as_tibble() %>%
    dplyr::rename(ensembl_gene_id_version = "geneID",
        ensembl_transcript_id_version = "txID") %>%
    left_join(biomart_gene, by = "ensembl_gene_id_version")

# Undo the name change
res_stageR <- res_stageR %>%
    mutate(ensembl_gene_id_version = sub('_', ':', ensembl_gene_id_version))

# Add proportions
props <- proportions(dmds) %>%
    tidyr::gather("sample_id", "proportion", -gene_id, -feature_id) %>%
    mutate(sample_id = substr(sample_id, 2, nchar(sample_id))) %>%
    left_join(sample_info %>% select(sample_id, condition_temp), by = "sample_id") %>%
    group_by(gene_id, feature_id, condition_temp) %>%
    summarise(proportion = mean(proportion)) %>%
    ungroup() %>%
    mutate(condition_temp = paste0("proportion_", condition_temp)) %>%
    tidyr::spread(condition_temp, proportion) %>%
    dplyr::rename(ensembl_gene_id_version = "gene_id",
        ensembl_transcript_id_version = "feature_id")

# export
log_info("Writing results table to disc...")
res_stageR %>%
    left_join(props, by = c("ensembl_gene_id_version", "ensembl_transcript_id_version")) %>%
    dplyr::select("ensembl_gene_id_version", "ensembl_transcript_id_version",
        "mgi_symbol", "description", "gene_biotype", everything()) %>%
    write_csv(snakemake@output[["res"]])

log_success("Done")
