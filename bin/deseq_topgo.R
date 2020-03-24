#!/usr/bin/Rscript --no-restore --no-environ --no-save

source(".Rprofile")
library("readr")
library("dplyr")
library("purrr")
library("topGO")

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

# reference for the GO analysis are all expressed genes
all_genes <- res$ensembl_gene_id_version

# function to generate genelist
get.geneList <- function(signif) {
    gl <- as.factor(as.integer(all_genes %in% signif$ensembl_gene_id_version))
    names(gl) <- all_genes
    # remove version
    names(gl) <- unlist(lapply(stringr::str_split(names(gl), "[.]"), "[[",1))
    gl
}

# function to actually create topGO object
make.topGO <- function(geneList, description) {
    new("topGOdata",
        description = description,
        ontology = "BP",
        allGenes = geneList,
        nodeSize = 10,
        annot = annFUN.org,
        ID = "ensembl",
        mapping = "org.Mm.eg")
}

# conveniance function to get results from topGO
get.results <- function(topGO) {
    topGOres <- runTest(topGO, algorithm = "parentchild", statistic = "fisher")
    GenTable(topGO, p = topGOres, orderBy = "p", ranksOf = "p", topNodes = 10) %>%
        as_tibble() %>%
        mutate(p = as.double(p))
}

# run the analysis
res_go <- res %>%
    filter(get(pvalue_name) < .05) %>%
    dplyr::select(ensembl_gene_id_version, log2FoldChange) %>%
    mutate(regulation = if_else(log2FoldChange > 0 , "up", "down")) %>%
    group_by(regulation) %>%
    tidyr::nest() %>%
    mutate(data = map(data, get.geneList)) %>%
    mutate(topGO = map2(data, as.character(regulation), make.topGO)) %>%
    mutate(res = map(topGO, get.results)) %>%
    mutate(n = map(data, function(x) sum(as.integer(as.character(x))))) %>%
    mutate(n = unlist(n))

# export
res_go$res %>%
    set_names(res_go$regulation) %>%
    saveRDS(snakemake@output[[1]])
