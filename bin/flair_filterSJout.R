#!/usr/bin/Rscript --no-restore --no-environ --no-save

# set libpaths to packrat local library
source(".Rprofile")
library("readr")
library("dplyr")

################################################################################
#
# author: Christoph Kiefer
# email: christophak@bmb.sdu.dk
#
################################################################################

df <- read_tsv(snakemake@input[[1]],
        col_names = FALSE,
        col_types = cols(
            X1 = col_character(),
            X2 = col_integer(),
            X3 = col_integer(),
            X4 = col_integer(),
            X5 = col_integer(),
            X6 = col_integer(),
            X7 = col_integer(),
            X8 = col_integer(),
            X9 = col_integer())) %>%
    filter(X7 >= snakemake@params[["threshold"]])

df %>%
    write_tsv(snakemake@output[[1]], col_names = FALSE)
