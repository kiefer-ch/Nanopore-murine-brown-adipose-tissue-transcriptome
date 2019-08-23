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
args = commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
    stop("At least one argument must be supplied (sample_info.csv).\n", call. = FALSE)
}

sample_df <- args[1]
HUB_URL <- args[2]

# expand variables
HUB <- "nanoporeibat_hub/hub/"
BW <- "nanoporeibat_hub/BW/"

# read sample info
sample_info <- read_csv(sample_df) %>%
    mutate_at(vars(matches("condition")), as.factor) %>%
    mutate(short_label = paste(sample_id, condition_temp)) %>%
    mutate(long_label = paste(sample_id, illumina, condition_temp)) %>%
    mutate(colour = brewer.pal(length(unique(condition_temp)), "Set1")[condition_temp]) %>%
    filter(!is.na(illumina))

# create directories
dir.create(file.path(HUB, "mm10"),
    showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(HUB, "../bw"),
    showWarnings = FALSE)

# create genomes.txt
cat("genome mm10\ntrackDb mm10/trackDb.txt",
    file = file.path(HUB, "genomes.txt"))

# create hub.txt
cat("hub nanoporeibat",
    "shortLabel illumina reads of murine iBAT",
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

    for (direction in c("fw", "rv")) {

        file <- sprintf("%s_%s.bw", row$illumina, direction)

        track <- c(
            sprintf("track %s_%s", row$sample_id, direction),
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

        file.copy(file.path("BW", file), file.path(HUB, "../bw", file))
    }
}

