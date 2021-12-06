
source(".Rprofile")
suppressPackageStartupMessages({
    library("logger")
    library("dplyr")
    library("dtplyr")
    library("data.table")
    library("purrr")
    library("GenomicAlignments")
})

################################################################################
#
# author: Christoph Kiefer
# email: christophak@bmb.sdu.dk
#
################################################################################

log_info(sprintf("Importing %s...", snakemake@input[[1]]))
bam <- scanBam(snakemake@input[[1]],
    param = ScanBamParam(
        what = c("rname"),
        flag = scanBamFlag(
            isSupplementaryAlignment = FALSE,
            isUnmappedQuery = FALSE)))[[1]] %>%
    lazy_dt()


log_info(sprintf("Processing %s...", snakemake@input[[1]]))
bam <- bam %>%
    mutate(rname = gsub("\\|.+", "", rname)) %>%
    group_by(rname) %>%
    summarise(NumReads = n())

log_info(sprintf("Writing %s to disk...", snakemake@output[[1]]))
bam %>%
    as_tibble() %>%
    dplyr::rename(Name = "rname") %>%
    saveRDS(snakemake@output[[1]])

log_success("Done.")
