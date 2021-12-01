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

get_map_type <- function(flag) {
    if_else(flag == 4, "unmapped",
            if_else(flag %in% c(0, 16), "primary",
                    if_else(flag %in% c(256, 275), "secondary", "supplementary")))
}


log_info("Importing data...")
df_tx <- c(snakemake@input[["teloprime_bam_tx"]],
           snakemake@input[["cdna_bam_tx"]],
           snakemake@input[["rna_bam_tx"]]) %>%
    set_names(tools::file_path_sans_ext(basename(.)))

names(df_tx) <- strsplit(names(df_tx), '_') %>%
    map(`[`, 1:2) %>%
    map(paste, collapse = '_')

log_info("Collapsing data...")
df_tx <- df_tx %>%
    map(read_rds) %>%
    map(mutate, type = get_map_type(flag)) %>%
    map(dplyr::select, -flag) %>%
    map(group_by, qname) %>%
    map(mutate, has_supplementary = any(type == "supplementary")) %>%
    map(ungroup) %>%
    bind_rows(.id = "sample") %>%
    tidyr::separate(sample, c("library", "barcode"), sep = '_') %>%
    mutate(library = factor(library, levels = c("illumina", "teloprime", "cdna", "rna")))

log_info("Writing to disc...")
df_tx %>%
    saveRDS(snakemake@output[[1]])

log_success("Done.")
