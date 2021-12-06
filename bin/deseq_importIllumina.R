
source(".Rprofile")
suppressPackageStartupMessages({
    library("logger")
    library("readr")
    library("dplyr")
    library("tximport")
    library("DESeq2")
})

################################################################################
#
# author: Christoph Kiefer
# email: christophak@bmb.sdu.dk
#
################################################################################

log_info("Generate tx2gene table...")
tx2g <- AnnotationDbi::loadDb(snakemake@input[["txdb"]])  %>%
    AnnotationDbi::select(.,
        keys =  AnnotationDbi::keys(., keytype = "GENEID"),
        keytype = "GENEID",
        columns = "TXNAME") %>%
    dplyr::select(TXNAME, GENEID)


log_info("Pprepare sample_info")
files <- tibble(path = snakemake@input$salmon_out) %>%
    mutate(illumina = path %>% dirname() %>% basename())

sample_info <- read_csv(snakemake@input[["sample_info"]],
        show_col_types = FALSE) %>%
    mutate_at(vars(matches("condition")), as.factor) %>%
    filter(!is.na(ont)) %>%
    select(sample_id, condition_temp, illumina) %>%
    left_join(files, "illumina")


log_info("Create and output dds...")
if(snakemake@params$type == "transcript") {
    tx_out <- TRUE
} else if(snakemake@params$type == "gene") {
    tx_out <- FALSE
}

tximport(files = sample_info$path,
        type = "salmon",
        tx2gene = tx2g,
        txOut = tx_out) %>%
    DESeqDataSetFromTximport(.,
        colData = tibble::column_to_rownames(as.data.frame(sample_info), "sample_id"),
        design = ~condition_temp) %>%
    saveRDS(snakemake@output[[1]])


log_success("Done")

