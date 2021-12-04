
source(".Rprofile")
suppressPackageStartupMessages({
    library("logger")
    library("dplyr")
    library("purrr")
    library("GenomicAlignments")
})

################################################################################
#
# author: Christoph Kiefer
# email: christophak@bmb.sdu.dk
#
# some code taken from https://github.com/csoneson/NativeRNAseqComplexTranscriptome/blob/master/Rscripts/get_nbr_reads.R
#
################################################################################

# function definitions
get_cigar <- function(cigar) {

    if(is.na(cigar)) {

        out <- list(NA_character_)

    } else {

        lengths = explodeCigarOpLengths(cigar)
        values = explodeCigarOps(cigar)

        out <- relist(paste0(unlist(lengths), unlist(values)), values)
    }

    out

}

get_opts <- function(cigar, opts) {
    # the next line is from https://github.com/csoneson/NativeRNAseqComplexTranscriptome/blob/master/Rscripts/get_nbr_reads.R
    sum(as.numeric(gsub(paste0(opts, "$"), "", cigar)), na.rm = TRUE)
}


log_info(sprintf("Importing %s...", snakemake@input[[1]]))
bam <- scanBam(snakemake@input[[1]],
    param = ScanBamParam(
        what = c("qname", "flag", "rname", "qwidth", "cigar", "strand", "pos")))[[1]]


log_info(sprintf("Processing %s...", snakemake@input[[1]]))
bam <- bam %>%
    as_tibble() %>%
    rowwise() %>%
    mutate(cigar = get_cigar(cigar)) %>%
    mutate(n_I = get_opts(cigar, "I"),
        n_M = get_opts(cigar, "M"),
        n_D = get_opts(cigar, "D"),
        n_N = get_opts(cigar, "N")) %>%
    ungroup() %>%
    mutate(aligned = as.integer(n_I + n_M),
           rwidth = as.integer(n_M + n_D + n_N)) %>%
    select(-cigar, -n_I, -n_M, -n_D, -n_N)


log_info(sprintf("Writing %s to disk...", snakemake@output[[1]]))
bam %>%
    saveRDS(snakemake@output[[1]])

log_success("Done.")
