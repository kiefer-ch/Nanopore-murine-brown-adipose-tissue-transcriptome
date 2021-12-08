
source(".Rprofile")
suppressPackageStartupMessages({
    library("logger")
    library("dplyr")
    library("dtplyr")
    library("data.table")
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
    read_rds() %>%
    select(ensembl_transcript_id_version, transcript_length)

df_tx <- c(snakemake@input$coverage) %>%
    set_names(tools::file_path_sans_ext(basename(.)) %>%
        strsplit(., '_') %>%
        map(`[`, 1:2) %>%
        map(paste, collapse = '_')) %>%
    map(read_rds)


log_info("Collapse data...")
df_tx <- df_tx %>%
    bind_rows(.id = "sample") %>%
    lazy_dt() %>%
    mutate(
        rname = sub("\\|.+", "", rname),
        type = case_when(
            flag %in% c(0L, 16L)    ~ "primary",
            flag %in% c(256L, 275L) ~ "secondary",
            TRUE                    ~ "supplementary")) %>%
    dplyr::select(-flag) %>%
    group_by(sample, qname) %>%
    mutate(has_supplementary = any(type == "supplementary")) %>%
    mutate(category = case_when(
        type %in% c("supplementary", "unmapped")    ~ type,
        has_supplementary                           ~ "primary_with_supplementary",
        TRUE                                        ~ "primary_wo_supplementary")) %>%
    ungroup() %>%
    select(-qname)
    left_join(biomart, by = c("rname" = "ensembl_transcript_id_version")) %>%
    mutate(coverage = coverage / transcript_length) %>%
    select(-has_supplementary) %>%
    as_tibble() %>%
    tidyr::separate(sample, c("library", "barcode"))


log_info("Writing to disc...")
df_tx %>%
    saveRDS(snakemake@output[[1]])


log_success("Done.")
