
source(".Rprofile")
suppressPackageStartupMessages({
    library("logger")
    library("dplyr")
    library("readr")
    library("purrr")
})

################################################################################
#
# author: Christoph Kiefer
# email: christophak@bmb.sdu.dk
#
################################################################################

log_info("Importing data...")
biomart <- snakemake@input[["biomaRt_tx"]] %>%
    read_rds()


df_tx <- c(snakemake@input$coverage) %>%
    set_names(tools::file_path_sans_ext(basename(.)) %>%
                  strsplit(., '_') %>%
                  map(`[`, 1:2) %>%
                  map(paste, collapse = '_')) %>%
    map(read_rds)


log_info("Collapse data...")
df_tx <- df_tx %>%
    map(mutate, type = case_when(
        flag %in% c(0L, 16L)    ~ "primary",
        flag %in% c(256L, 275L) ~ "secondary",
        TRUE                    ~ "supplementary")) %>%
    map(filter, type == "primary") %>%
    map(dplyr::select, -qname, -flag) %>%
    map(tidyr::separate, col = seqnames, into = "ensembl_transcript_id_version",
        sep = "\\|", extra = "drop") %>%
    map(left_join, y = biomart, by = "ensembl_transcript_id_version") %>%
    map(mutate, coverage = coverage / transcript_length) %>%
    map(dplyr::select, coverage, transcript_length, type, qwidth) %>%
    bind_rows(.id = "sample") %>%
    tidyr::separate(sample, c("library", "barcode"))


log_info("Writing to disc...")
df %>%
    saveRDS(snakemake@output[[1]])


log_success("Done.")
