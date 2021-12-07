
source(".Rprofile")
suppressPackageStartupMessages({
    library("logger")
    library("dplyr")
    library("dtplyr")
    library("data.table")
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
bam <- scanBam(snakemake@input[[1]],
        param = ScanBamParam(
           what = c("qname", "flag", "rname", "pos")))[[1]] %>%
    lazy_dt()

readGAlignments(snakemake@input[[1]],
    use.names = TRUE,
    param = ScanBamParam(tag = c("NM"),
        scanBamFlag(isUnmappedQuery = NA),
        what = c("qname","flag", "rname", "pos"))) %>%
    as_tibble() %>%
    select(qname, flag, qwidth, cigar, seqnames) %>%
    mutate(lengths = explodeCigarOpLengths(cigar),
        values = explodeCigarOps(cigar)) %>%
    mutate(cigar = relist(paste0(unlist(lengths), unlist(values)), values)) %>%
    select(-lengths, -values) %>%
    rowwise() %>%
    mutate(n_D = get_opts(cigar, "D"),
        n_M = get_opts(cigar, "M")) %>%
    ungroup() %>%
    mutate(coverage = n_D + n_M) %>%
    select(-cigar, -n_D, -n_M) %>%
    saveRDS(snakemake@output[[1]])
