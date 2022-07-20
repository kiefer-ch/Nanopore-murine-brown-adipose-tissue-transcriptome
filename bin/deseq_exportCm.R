
source(".Rprofile")
suppressPackageStartupMessages({
    library("readr")
    library("dplyr")
    library("DESeq2")
    library("logger")
})

################################################################################
#
# author: Christoph Kiefer
# email: christophak@bmb.sdu.dk
#
################################################################################

# import data
log_info("Import data...")
dds <- readRDS(snakemake@input[["dds"]]) %>%
     estimateSizeFactors()


log_info("Exporting raw counts...")
counts(dds, normalized = FALSE) %>%
    as_tibble(rownames = "id") %>%
    write_csv(snakemake@output[["cts"]])


# simple median gene normalisation
log_info("Exporting normalised data...")
counts(dds, normalized = TRUE) %>%
    as_tibble(rownames = "id") %>%
    write_csv(snakemake@output[["ntd"]])


# tpm
if(snakemake@wildcards$dataset == "illumina") {
    # https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/
    fpkmToTpm <- function(fpkm) {
        exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
    }

    log_info("Exporting TPM...")
    tpm <- fpkm(dds, robust = TRUE)
    tpm <- apply(tpm, 2, fpkmToTpm)

    tpm %>%
        as_tibble(rownames = "id") %>%
        write_csv(snakemake@output[["tpm"]])
}


log_success("Done.")
