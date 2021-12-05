
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
df_tx <- c(snakemake@input[[1]]) %>%
    set_names(tools::file_path_sans_ext(basename(.)) %>%
        strsplit(., '_') %>%
        map(`[`, 1:3) %>%
        map(paste, collapse = '_')) %>%
    map(read_rds)

log_info("Collapsing data...")
df_tx <- df_tx %>%
    map(mutate, type = case_when(
        flag == 4L              ~ "unmapped",
        flag %in% c(0L, 16L)    ~ "primary",
        flag %in% c(256L, 275L) ~ "secondary",
        TRUE                    ~ "supplementary")) %>%
    map(dplyr::select, -flag) %>%
    map(group_by, qname) %>%
    map(mutate, has_supplementary = any(type == "supplementary")) %>%
    map(mutate, category = case_when(
        type %in% c("supplementary", "unmapped")    ~ type,
        has_supplementary                           ~ "primary_with_supplementary",
        TRUE                                        ~ "primary_wo_supplementary")) %>%
    map(ungroup) %>%
    bind_rows(.id = "sample") %>%
    select(-has_supplementary) %>%
    tidyr::separate(sample, c("library", "type", "barcode"), sep = '_') %>%
    select(-type)

log_info("Writing to disc...")
df_tx %>%
    saveRDS(snakemake@output[[1]])

log_success("Done.")
