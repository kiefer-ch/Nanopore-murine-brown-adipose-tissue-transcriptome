
source(".Rprofile")
suppressPackageStartupMessages({
    library("logger")
    library("dplyr")
    library("purrr")
    library("readr")
})

################################################################################
#
# author: Christoph Kiefer
# email: christophak@bmb.sdu.dk
#
################################################################################

log_info("Pprepare sample_info...")
if(snakemake@wildcards$dataset %in% c("teloprime", "cdna")) {

    if(snakemake@wildcards$dataset == "teloprime") {
        handle <- "ont"
    } else {
        handle <- "cdna"
    }

    sample_info <- read_csv(snakemake@input[["sample_info"]],
                            show_col_types = FALSE) %>%
        mutate_at(vars(matches("condition")), as.factor) %>%
        filter(!is.na(ont)) %>%
        dplyr::select(sample_id, condition_temp, all_of(handle))

} else if (snakemake@wildcards$dataset == "rna") {

    handle <- "rna"

    sample_info <- tibble(sample_id = c("warm", "cold"),
        condition_temp = as.factor(c(22, 4)),
        rna = c("rt", "cool"))
}


log_info("Prepare counts...")
counts <- snakemake@input$counts %>%
    set_names(basename(.) %>%
                  sub("_quant.tsv", "", .)) %>%
    map(read_tsv, col_types = "ci") %>%
    bind_rows(.id = "id") %>%
    tidyr::separate(id, c("trash", "trash2", "barcode")) %>%
    left_join(sample_info %>% dplyr::select(sample_id, handle), by = c("barcode" = handle)) %>%
    dplyr::select(-trash, -trash2, -barcode)



if (snakemake@params$type == "gene") {

    log_info("Generate tx2gene table...")
    tx2g <- AnnotationDbi::loadDb(snakemake@input[["txdb"]])  %>%
        AnnotationDbi::select(.,
                              keys =  AnnotationDbi::keys(., keytype = "GENEID"),
                              keytype = "GENEID",
                              columns = "TXNAME") %>%
        dplyr::select(TXNAME, GENEID)

    log_info("Aggregating genes...")
    counts <- counts %>%
        left_join(tx2g, by = c("Name" = "TXNAME")) %>%
        dplyr::select(-Name) %>%
        dplyr::rename(Name = "GENEID") %>%
        group_by(sample_id, Name) %>%
        summarise_if(is.numeric, sum)
}


counts <- counts %>%
    tidyr::pivot_wider(names_from = sample_id,
                       values_from = NumReads,
                       values_fill = 0) %>%
    dplyr::select(Name, sample_info$sample_id) %>%
    as.data.frame() %>%
    tibble::column_to_rownames("Name") %>%
    data.matrix()


log_info("Create dds...")
dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts,
    colData = tibble::column_to_rownames(as.data.frame(sample_info), "sample_id"),
    design = ~condition_temp)


log_info("Writing to disc...")
saveRDS(dds, snakemake@output[[1]])


log_success("Done")
