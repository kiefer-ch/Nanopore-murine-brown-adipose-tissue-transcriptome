
source(".Rprofile")
suppressPackageStartupMessages({
    library("dplyr")
    library("GenomicAlignments")
})

################################################################################
#
# author: Christoph Kiefer
# email: christophak@bmb.sdu.dk
#
################################################################################

get_opts <- function(cigar, opts) {
    # the next line is from https://github.com/csoneson/NativeRNAseqComplexTranscriptome/blob/master/Rscripts/get_nbr_reads.R
    as.integer(sum(suppressWarnings(as.numeric(gsub(paste0(opts, "$"), "", cigar))), na.rm = TRUE))
}

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
