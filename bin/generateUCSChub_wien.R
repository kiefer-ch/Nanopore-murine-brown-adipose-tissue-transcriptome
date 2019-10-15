#!/usr/bin/Rscript --no-restore --no-environ --no-save

rm(list = ls())

###############################################################################
#
# author: Christoph Kiefer
# email: christophak@bmb.sdu.dk
#
################################################################################

library("dplyr")
library("readr")
library("RColorBrewer")

################################################################################
# read command line args
sample_df <- "../sample_info/sampleInfo.csv"
HUB_URL <- "http://bioinformatik.sdu.dk/solexa/webshare/christoph/nanoporeibat_ont_hub"

# expand variables
HUB <- "nanoporeibat_ont_hub/hub/"
BW <- "nanoporeibat_ont_hub/BW/"

# read sample info
sample_info <- read_csv(sample_df) %>%
    mutate_at(vars(matches("condition")), as.factor) %>%
    mutate(short_label = paste(sample_id, condition_temp)) %>%
    mutate(long_label = paste(sample_id, illumina, condition_temp)) %>%
    mutate(colour = brewer.pal(length(unique(condition_temp)), "Set1")[condition_temp]) %>%
    filter(!is.na(ont)) %>%
    arrange(condition_temp)

# create directories
dir.create(file.path(HUB, "mm10"),
    showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(HUB, "../bw"),
    showWarnings = FALSE)

# create genomes.txt
cat("genome mm10\ntrackDb mm10/trackDb.txt",
    file = file.path(HUB, "genomes.txt"))

# create hub.txt
cat("hub nanoporeibat_ont",
    "shortLabel ont reads of murine iBAT",
    "longLabel RNA from BAT of cold treated and control mice.",
    "genomesFile genomes.txt",
    "email  christophak@bmb.sdu.dk",
    "descriptionUrl",
    sep = '\n',
    file = file.path(HUB, "hub.txt"))

# create trackDb.txt
suppressWarnings(file.remove(file.path(HUB, "hub", "mm10", "trackDb.txt")))
for (i in 1:nrow(sample_info)) {

    row <- sample_info[i,]

    file <- sprintf("%s.bw", row$ont)

    track <- c(
        sprintf("track %s", row$sample_id),
        sprintf("bigDataUrl %s", paste0(HUB_URL, "/bw/", file)),
        sprintf("shortLabel %s", row$short_label),
        sprintf("longLabel %s", row$long_label),
        "type bigWig",
        "visibility full",
        "autoScale on",
        sprintf("color %s", col2rgb(row$colour) %>%
            as.vector() %>%
            paste(collapse = ',')),
        "")

    cat(track,
        file = file.path(HUB, "mm10", "trackDb.txt"),
        sep = '\n',
        append = TRUE)

    file.copy(file.path("../BW/bw_ont", file), file.path(HUB, "../bw", file))
}

