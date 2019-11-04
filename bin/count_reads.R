#!/usr/bin/Rscript --no-restore --no-environ --no-save

# set libpaths to packrat local library
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
message("Preparing genome bam...")
bam <- readGAlignments(snakemake@input[["bam_genome"]],
    use.names = TRUE,
    param = ScanBamParam(tag = c("NM"),
        what = c("qname","flag", "rname", "pos", "mapq")))

get_opts <- function(cigar, opts) {
    sum(as.numeric(gsub(paste0(opts, "$"), "", cigar)), na.rm = TRUE)
}

df_genome <- bam %>%
    as_tibble() %>%
    select(qname, flag, qwidth, cigar) %>%
    mutate(lengths = explodeCigarOpLengths(cigar),
        values = explodeCigarOps(cigar)) %>%
    mutate(cigar = relist(paste0(unlist(lengths), unlist(values)), values)) %>%
    select(-lengths, -values) %>%
    rowwise() %>%
    mutate(n_I = get_opts(cigar, "I"),
        n_M = get_opts(cigar, "M")) %>%
    mutate(aligned = n_I + n_M) %>%
    select(-cigar, -n_I, -n_M) %>%
    group_by(qname) %>%
    summarise(n_primary_genome = sum(flag %in% c(0, 16)),
        n_supplementary_genome = sum(flag %in% c(2048, 2064)),
        qwidth_genome = max(qwidth),
        qaligned_genome = max(aligned))

# transcriptome
message("Preparing transcriptome bam...")
bam <- readGAlignments(snakemake@input[["bam_tx"]],
    use.names = TRUE,
    param = ScanBamParam(tag = c("NM"),
        what = c("qname","flag", "rname", "pos", "mapq")))

df_transcriptome <- bam %>%
    as_tibble() %>%
    select(qname, flag, qwidth, cigar) %>%
    mutate(lengths = explodeCigarOpLengths(cigar),
        values = explodeCigarOps(cigar)) %>%
    mutate(cigar = relist(paste0(unlist(lengths), unlist(values)), values)) %>%
    select(-lengths, -values) %>%
    rowwise() %>%
    mutate(n_I = get_opts(cigar, "I"),
        n_M = get_opts(cigar, "M")) %>%
    mutate(aligned = n_I + n_M) %>%
    select(-cigar, -n_I, -n_M) %>%
    group_by(qname) %>%
    summarise(n_primary_tx = sum(flag %in% c(0, 16)),
        n_secondary_tx = sum(flag %in% c(256, 272)),
        n_supplementary_tx = sum(flag %in% c(2048, 2064)),
        qwidth_tx = max(qwidth),
        qaligned_tx = max(aligned))

# merge and export
message("Exporting...")
df_genome %>%
    full_join(df_transcriptome, by = "qname") %>%
    saveRDS(snakemake@output[[1]])

