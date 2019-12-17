#!/usr/bin/Rscript --no-restore --no-environ --no-save

# set libpaths to packrat local library

save.image("coverage.RData")

message("Load packages...")
source("packrat/init.R")
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("purrr"))
suppressPackageStartupMessages(library("GenomicAlignments"))

################################################################################
#
# author: Christoph Kiefer
# email: christophak@bmb.sdu.dk
#
# some code taken from https://github.com/csoneson/NativeRNAseqComplexTranscriptome/blob/master/Rscripts/get_nbr_reads.R
#
################################################################################

# genome
message(sprintf("Loading %s...", snakemake@input[[1]]))
bam <- readGAlignments(snakemake@input[[1]],
    use.names = TRUE,
    param = ScanBamParam(tag = c("NM"),
        scanBamFlag(isUnmappedQuery = NA),
        what = c("qname","flag", "rname", "pos")))

get_opts <- function(cigar, opts) {
    sum(as.numeric(gsub(paste0(opts, "$"), "", cigar)), na.rm = TRUE)
}

message(sprintf("Processing %s...", snakemake@input[[1]]))
bam %>%
    as_tibble() %>%
    select(qname, flag, qwidth, cigar, seqnames) %>%
    mutate(lengths = explodeCigarOpLengths(cigar),
        values = explodeCigarOps(cigar)) %>%
    mutate(cigar = relist(paste0(unlist(lengths), unlist(values)), values)) %>%
    select(-lengths, -values) %>%
    rowwise() %>%
    mutate(n_I = get_opts(cigar, "I"),
        n_M = get_opts(cigar, "M")) %>%
    ungroup() %>%
    mutate(aligned = n_I + n_M) %>%
    select(-cigar, -n_I, -n_M) %>%
    saveRDS(snakemake@output[[1]])
