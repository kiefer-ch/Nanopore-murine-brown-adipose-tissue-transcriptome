
source(".Rprofile")
suppressPackageStartupMessages({
    library("logger")
    library("readr")
    library("dplyr")
    library("tximport")
    library("DRIMSeq")
})

################################################################################
#
# author: Christoph Kiefer
# email: christophak@bmb.sdu.dk
#
################################################################################

log_info("Pprepare sample_info")
files <- tibble(path = snakemake@input$salmon_out) %>%
    mutate(illumina = path %>% dirname() %>% basename())

sample_info <- read_csv(snakemake@input[["sample_info"]],
        show_col_types = FALSE) %>%
    mutate_at(vars(matches("condition")), as.factor) %>%
    filter(!is.na(ont)) %>%
    dplyr::select(sample_id, condition_temp, illumina) %>%
    left_join(files, "illumina")



log_info("Generate tx2gene table...")
tx2g <- AnnotationDbi::loadDb(snakemake@input[["txdb"]])  %>%
    AnnotationDbi::select(.,
        keys =  AnnotationDbi::keys(., keytype = "GENEID"),
        keytype = "GENEID",
        columns = "TXNAME") %>%
    dplyr::select(TXNAME, GENEID)


if (snakemake@wildcards$method == "flair") {

   tx2g$TXNAME <-  paste(tx2g$TXNAME, tx2g$GENEID, sep = "_")

}

log_info("Import reads...")
txi <- tximport(files = sample_info$path,
    type = "salmon",
    tx2gene = tx2g,
    txOut = TRUE,
    countsFromAbundance = "dtuScaledTPM")

cts <- txi$counts
cts <- cts[rowSums(cts) > 0, ]
colnames(cts) <- sample_info$sample_id


log_info("Create txdf...")
txdf <- tx2g %>%
    group_by(GENEID) %>%
    mutate(ntx = row_number()) %>%
    ungroup()

txdf <- txdf[match(rownames(cts), txdf$TXNAME),]

counts <- cts %>%
    as_tibble(rownames = "TXNAME") %>%
    left_join(txdf, by = "TXNAME") %>%
    dplyr::rename(gene_id = "GENEID",
        feature_id = "TXNAME") %>%
    as.data.frame()


log_info("Create dnds...")
dmds <- dmDSdata(counts = counts, samples = as.data.frame(sample_info))


log_info("Write to disc...")
saveRDS(cts, snakemake@output[[1]])


log_success("Done.")
