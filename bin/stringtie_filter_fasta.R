#!/usr/bin/Rscript --no-restore --no-environ --no-save

source("packrat/init.R")
library("readr")
library("dplyr")

################################################################################
#
# author: Christoph Kiefer
# email: christophak@bmb.sdu.dk
#
# The read_gtf function is from Paul Klemm
#
################################################################################

read_gtf <- function(file) {
    read_delim(file = file,
        delim = "\t",
        comment = "#",
        na = c('.'),
        col_names = c("sequence", "source", "feature", "start", "end", "score",
            "strand", "phase", "attributes"),
        col_types = cols(
            sequence = col_character(),
            source = col_character(),
            feature = col_character(),
            start = col_integer(),
            end = col_integer(),
            score = col_character(),
            strand = col_character(),
            phase = col_character(),
            attributes = col_character()),
        progress = FALSE)
}

write_gtf <- function(gtf, path) {
    withr::with_options(c(scipen = 10),
        write.table(gtf, path,
            col.names = FALSE, row.names = FALSE,
            na = '.', quote = FALSE, sep = '\t'))
}

# read gtf
gtf <- read_gtf(snakemake@input[[1]])

# filter NA strand
gtf <- gtf %>%
    filter(!is.na(strand))

# write filtered gtf
write_gtf(gtf, snakemake@output[[1]])
