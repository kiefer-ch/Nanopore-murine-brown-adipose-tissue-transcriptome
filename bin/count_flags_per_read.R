#!/usr/bin/Rscript --no-restore --no-environ --no-save

# set libpaths to packrat local library
source("packrat/init.R")

library("dplyr")
library("GenomicAlignments")

################################################################################
#
# author: Christoph Kiefer
# email: christophak@bmb.sdu.dk
#
################################################################################

# genome
bam <- readGAlignments("../BAM/bam_ont/barcode01.bam",
    use.names = TRUE,
    param = ScanBamParam(tag = c("NM"),
        what = c("qname","flag", "rname", "pos", "mapq")))

bam <- bam %>%
    as_tibble()

bam %>%
    `[`(1:10000,) %>%
    select(qname, flag) %>%
    group_by(qname) %>%
    summarise(has_primary_genome = any(flag %in% c(0, 16)),
        has_supplementary_genome = any(flag %in% c(2048, 2064)))
