
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

# function definitions
get_cigar <- function(cigar) {
    if(is.na(cigar)) {
        out <- list(NA_character_)
    } else {
        lengths = explodeCigarOpLengths(cigar)
        values = explodeCigarOps(cigar)
        out <- unlist(relist(paste0(unlist(lengths), unlist(values)), values))
    }
    out
}

get_opts <- function(cigar, opts) {
    # the next line is from https://github.com/csoneson/NativeRNAseqComplexTranscriptome/blob/master/Rscripts/get_nbr_reads.R
    as.integer(sum(suppressWarnings(as.numeric(gsub(paste0(opts, "$"), "", cigar))), na.rm = TRUE))
}


log_info(sprintf("Importing %s...", snakemake@input[[1]]))
# skip unmapped reads
bam <- scanBam(snakemake@input[[1]],
        param = ScanBamParam(
            scanBamFlag(isUnmappedQuery = FALSE),
            what = c("flag", "rname", "cigar")))[[1]] %>%
    lazy_dt()


log_info(sprintf("Processing %s...", snakemake@input[[1]]))
bam <- bam %>%
    mutate(cigar = map(cigar, get_cigar)) %>%
    mutate(n_M = map_int(cigar, get_opts, opts = "M"),
           n_D = map_int(cigar, get_opts, opts = "D")) %>%
    mutate(coverage = n_D + n_M) %>%
    select(-cigar, -n_M)


log_info(sprintf("Writing %s to disk...", snakemake@output[[1]]))
bam %>%
    as_tibble() %>%
    saveRDS(snakemake@output[[1]])

log_success("Done.")
