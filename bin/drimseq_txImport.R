
save.image()

source(".Rprofile")
suppressPackageStartupMessages({
    library("logger")
    library("readr")
    library("dplyr")
    library("tximport")
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

# generate tx2gene table
tx2g <- AnnotationDbi::loadDb(snakemake@input[["txdb"]])  %>%
    AnnotationDbi::select(.,
        keys =  AnnotationDbi::keys(., keytype = "GENEID"),
        keytype = "GENEID",
        columns = "TXNAME") %>%
    dplyr::select(TXNAME, GENEID)

# create and output dds
txi <- tximport(files = sample_info$path,
    type = "salmon",
    tx2gene = tx2g,
    txOut = TRUE,
    countsFromAbundance = "dtuScaledTPM")

cts <- txi$counts
cts <- cts[rowSums(cts) > 0, ]
colnames(cts) <- sample_info$sample_id

saveRDS(cts, snakemake@output[[1]])
