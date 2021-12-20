
source(".Rprofile")
suppressPackageStartupMessages({
    library("dplyr")
    library("logger")
})

################################################################################
#
# author: Christoph Kiefer
# email: christophak@bmb.sdu.dk
#
################################################################################
log_info("MakingTxDb from gtf...")
if (snakemake@params[["species"]] == "human") {
    species <- "Homo sapiens"
} else if (snakemake@params[["species"]] == "mouse") {
    species <- "Mus musculus"
} else {
    species <- NA_character_
}

GenomicFeatures::makeTxDbFromGFF(snakemake@input[[1]],
        format = "gtf",
        organism = species,
        circ_seqs = character()) %>%
    AnnotationDbi::saveDb(snakemake@output[[1]])

log_success("Done.")

