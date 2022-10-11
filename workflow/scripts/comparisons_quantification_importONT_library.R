
source(".Rprofile")
suppressPackageStartupMessages({
    library("logger")
    library("dplyr")
    library("purrr")
    library("readr")
    library("DESeq2")
})

################################################################################
#
# author: Christoph Kiefer
# email: christophak@bmb.sdu.dk
#
################################################################################

log_info("Pprepare sample_info...")
sample_info <- read_csv(snakemake@input[["sample_info"]],
                        show_col_types = FALSE) %>%
    mutate_at(vars(matches("condition")), as.factor) %>%
    filter(!is.na(ont)) %>%
    dplyr::select(ont, cdna, condition_temp) %>%
    tidyr::pivot_longer(cols = c(ont, cdna), names_to = "library",
        values_to = "barcode") %>%
    mutate(library = if_else(library == "ont", "teloprime", library)) %>%
    mutate(x1 = "flowcell1", x2 = "flowcell2") %>%
    tidyr::pivot_longer(cols = c(x1, x2), values_to = "flowcell") %>%
    select(-name) %>%
    tidyr::unite("sample_id", library, flowcell, barcode, remove = FALSE) %>%
    arrange(flowcell)


log_info("Pprepare count data...")
cdna <- snakemake@input$cdna %>%
    set_names(basename(.) %>%
                  sub("_quant.tsv", "", .)) %>%
    map(read_tsv, col_types = "ci") %>%
    bind_rows(.id = "id") %>%
    tidyr::pivot_wider(names_from = id, values_from = NumReads, values_fill = 0) %>%
    as.data.frame() %>%
    tibble::column_to_rownames("Name") %>%
    data.matrix() %>%
    DESeq2::DESeqDataSetFromMatrix(countData = .,
        colData = tibble::column_to_rownames(as.data.frame(sample_info %>% filter(library == "cdna")), "sample_id"),
        design = ~condition_temp * flowcell) %>%
    estimateSizeFactors() %>%
    counts(normalized = TRUE) %>%
    as_tibble(rownames = "id")

teloprime <- snakemake@input$teloprime %>%
    set_names(basename(.) %>%
                  sub("_quant.tsv", "", .)) %>%
    map(read_tsv, col_types = "ci") %>%
    bind_rows(.id = "id") %>%
    tidyr::pivot_wider(names_from = id, values_from = NumReads, values_fill = 0) %>%
    as.data.frame() %>%
    tibble::column_to_rownames("Name") %>%
    data.matrix() %>%
    DESeq2::DESeqDataSetFromMatrix(countData = .,
                                   colData = tibble::column_to_rownames(as.data.frame(sample_info %>% filter(library == "teloprime")), "sample_id"),
                                   design = ~condition_temp * flowcell) %>%
    estimateSizeFactors() %>%
    counts(normalized = TRUE) %>%
    as_tibble(rownames = "id")

df <- full_join(cdna, teloprime, by = "id")


log_info("Writing to disc...")
write_tsv(df, snakemake@output[[1]])


log_success("Done")
