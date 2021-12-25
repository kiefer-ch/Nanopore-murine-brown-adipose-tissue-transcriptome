
source(".Rprofile")
suppressPackageStartupMessages({
    library("logger")
    library("dplyr")
    library("readr")
    library("DRIMSeq")
    library("stageR")
})

################################################################################
#
# author: Christoph Kiefer
# email: christophak@bmb.sdu.dk
#
################################################################################

log_info("Importing data...")
dmds <- read_rds(snakemake@input$dmds)

biomart_gene <- read_rds(snakemake@input$biomart)

sample_info <- read_csv(snakemake@input[["sample_info"]],
        show_col_types = FALSE) %>%
    mutate_at(vars(matches("condition")), as.factor) %>%
    filter(!is.na(ont)) %>%
    dplyr::select(condition_temp, sample_id, cdna)


if (snakemake@wildcards$dataset == "cdna") {
    sample_info <- sample_info %>%
        select(condition_temp, sample_id = "cdna")
} else {
    sample_info <- sample_info %>%
        select(-cdna)
}


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
    dplyr::rename(gene_id = "geneID",
        transcript_id = "txID")


# Add proportions
props <- proportions(dmds) %>%
    tidyr::gather("sample_id", "proportion", -gene_id, -feature_id)

if (snakemake@wildcards$dataset == "illumina") {
    props <- props %>%
    mutate(sample_id = substr(sample_id, 2, nchar(sample_id)))
}

props <- props %>%
    left_join(sample_info, by = "sample_id") %>%
    group_by(gene_id, feature_id, condition_temp) %>%
    summarise(proportion = mean(proportion)) %>%
    ungroup() %>%
    mutate(condition_temp = paste0("proportion_", condition_temp)) %>%
    tidyr::spread(condition_temp, proportion) %>%
    dplyr::rename(ensembl_gene_id_version = "gene_id",
        ensembl_transcript_id_version = "feature_id")

res_stageR <- res_stageR %>%
    left_join(props, by = c("gene_id" = "ensembl_gene_id_version",
                        "transcript_id" = "ensembl_transcript_id_version"))


if (!is.null(snakemake@input$tmap)) {
    tmap <- read_tsv(snakemake@input$tmap, col_types = "ccfccidddici", na = "-") %>%
        select(ref_gene_id, ref_id, qry_gene_id, qry_id)

    res_stageR <- res_stageR %>%
        left_join(tmap,
                  by = c("gene_id" = "qry_gene_id", "transcript_id" = "qry_id")) %>%
        mutate(gene_id = if_else(is.na(ref_gene_id), gene_id, ref_gene_id),
               transcript_id = if_else(is.na(ref_id), transcript_id, ref_id)) %>%
        select(-ref_gene_id, -ref_id)
}

res_stageR <- res_stageR %>%
    left_join(biomart_gene, by = c("gene_id" = "ensembl_gene_id_version"))


# Undo the name change
res_stageR <- res_stageR %>%
    mutate(gene_id = sub('_', ':', gene_id))


log_info("Writing results table to disc...")
res_stageR %>%
    dplyr::select("gene_id", "transcript_id",
        "mgi_symbol", "description", "gene_biotype", everything()) %>%
    write_csv(snakemake@output[[1]])

log_success("Done")
