#!/usr/bin/Rscript --no-restore --no-environ --no-save

# set libpaths to packrat local library
source(".Rprofile")

suppressPackageStartupMessages(library("dplyr"))

################################################################################
#
# author: Christoph Kiefer
# email: christophak@bmb.sdu.dk
#
################################################################################
GenomicFeatures::makeTxDbFromGFF(snakemake@input[[1]],
        format = "gtf", circ_seqs = character()) %>%
    AnnotationDbi::saveDb(snakemake@output[[1]])
