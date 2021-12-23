
source(".Rprofile")
suppressPackageStartupMessages({
    library("logger")
    library("readr")
    library("dplyr")
    library("DRIMSeq")
})

################################################################################
#
# author: Christoph Kiefer
# email: christophak@bmb.sdu.dk
#
################################################################################

log_info("Create txdf...")
txdf <- AnnotationDbi::loadDb(snakemake@input[["txdb"]])  %>%
    AnnotationDbi::select(.,
        keys =  AnnotationDbi::keys(., keytype = "GENEID"),
        keytype = "GENEID",
        columns = "TXNAME") %>%
    dplyr::select(TXNAME, GENEID) %>%
    group_by(GENEID) %>%
    mutate(ntx = row_number()) %>%
    ungroup()


log_info("Import counts...")
cts <- snakemake@input[["counts"]] %>%
    purrr::set_names(basename(.) %>%
        sub("_quant.tsv", "", .) %>%
            sub("^cdna_merged_", "", .)) %>%
    purrr::map(read_tsv, col_types = "ci") %>%
    bind_rows(.id = "barcode") %>%
    dplyr::rename(TXNAME = "Name") %>%
    tidyr::pivot_wider(names_from = barcode, values_from = "NumReads", values_fill = 0)


txdf <- txdf[match(cts$TXNAME, txdf$TXNAME),]


counts <- cts %>%
    left_join(txdf, by = "TXNAME") %>%
    dplyr::rename(gene_id = "GENEID",
        feature_id = "TXNAME") %>%
    as.data.frame()


log_info("Import sample info...")
sample_info <- read_csv(snakemake@input[["sample_info"]],
                        show_col_types = FALSE) %>%
    mutate_at(vars(matches("condition")), as.factor) %>%
    filter(!is.na(ont)) %>%
    dplyr::select(condition_temp, sample_id = "cdna")


log_info("Create dmds")
dmds <- dmDSdata(counts = counts, samples = as.data.frame(sample_info))



log_info("Write to disc...")
saveRDS(cts, snakemake@output[[1]])


log_success("Done.")
