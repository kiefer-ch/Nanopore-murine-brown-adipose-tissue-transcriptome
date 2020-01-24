#!/usr/bin/Rscript --no-restore --no-environ --no-save

source("packrat/init.R")
library("readr")
library("dplyr")
library("ReactomePA")

################################################################################
#
# author: Christoph Kiefer
# email: christophak@bmb.sdu.dk
#
################################################################################

# tx and gene level
if (snakemake@params[["level"]] == "txlevel") {
    id <- "ensembl_transcript_id_version"
} else if (snakemake@params[["level"]] == "genelevel") {
    id <- "ensembl_gene_id_version"
}

# shrunken or not
if (snakemake@params[["shrink"]] == "apeglm") {
    pvalue_name <- "svalue"
} else if (snakemake@params[["shrink"]] == "noShrink") {
    pvalue_name <- "padj"
}

# import data
res <- read_csv(snakemake@input[[1]])
ensembl = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# annotate with entrez ids
df_reactome <- res %>%
    pull(id) %>%
    biomaRt::getBM(
        attributes = c("entrezgene_id", id),
        filters = id,
        values = .,
        mart = ensembl) %>%
    left_join(res, by = id) %>%
    tidyr::drop_na(entrezgene_id) %>%
    mutate(entrezgene_id = as.character(entrezgene_id))

# reference are all expressed genes
geneList_reference <- df_reactome %>%
    pull(entrezgene_id)

# named vector with lfc and entrez ids
# significant: 5%
geneList_signif <- df_reactome %>%
    filter(get(pvalue_name) < .05) %>%
    dplyr::select(log2FoldChange, entrezgene_id) %>%
    tibble::deframe()
saveRDS(geneList_signif, snakemake@output[[2]])

# perform pathway analysis
res_reactome <- enrichPathway(geneList_signif,
        organism = "mouse",
        readable = TRUE,
        universe = geneList_reference)

# save to disc
saveRDS(res_reactome, snakemake@output[[1]])
