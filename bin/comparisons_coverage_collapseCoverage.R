
source(".Rprofile")
suppressPackageStartupMessages(library("logger"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("purrr"))

################################################################################
#
# author: Christoph Kiefer
# email: christophak@bmb.sdu.dk
#
################################################################################

log_info("Importing data...")
biomart <- snakemake@input[["biomaRt_tx"]] %>%
    read_rds()


get_map_type <- function(flag) {
    if_else(flag %in% c(0, 16), "primary",
            if_else(flag %in% c(256, 275), "secondary", "supplementary"))
}

log_info("Collapse data...")
prepare_dataset2 <- function(files) {
    files %>%
        as.list() %>%
        set_names(basename(files)) %>%
        map(read_rds) %>%
        map(mutate, type = get_map_type(flag)) %>%
        map(dplyr::select, -qname, -flag, -qwidth) %>%
        map(tidyr::separate, col = seqnames, into = "ensembl_transcript_id_version",
            sep = "\\|", extra = "drop") %>%
        map(left_join, y = biomart, by = "ensembl_transcript_id_version") %>%
        bind_rows(.id = "sample") %>%
        mutate(coverage = coverage / transcript_length)
}


df <- list(snakemake@input[["coverage_teloprime"]],
           snakemake@input[["coverage_cdna"]],
           snakemake@input[["coverage_rna"]]) %>%
    map(prepare_dataset2) %>%
    map(dplyr::select, coverage, transcript_length, type, ensembl_transcript_id_version) %>%
    map(filter, type == "primary") %>%
    set_names(c("teloprime", "cdna", "rna")) %>%
    bind_rows(.id = "dataset") %>%
    tidyr::drop_na()


log_info("Writing to disc...")
df %>%
    saveRDS(snakemake@output[[1]])


log_success("Done.")

