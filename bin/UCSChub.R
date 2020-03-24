#!/usr/bin/Rscript --no-restore --no-environ --no-save

###############################################################################
#
# author: Christoph Kiefer
# email: christophak@bmb.sdu.dk
#
################################################################################

source(".Rprofile")
library("dplyr")
library("readr")
library("RColorBrewer")

################################################################################

# expand variables
hub <- "nanoporeibat_hub/hub/"
bw <- "nanoporeibat_hub/bw/"

# read sample info
sample_info <- read_csv(snakemake@input[["sample_info"]]) %>%
    filter(!is.na(ont)) %>%
    select(sample_id, condition_temp) %>%
    mutate(condition_temp = as.factor(condition_temp)) %>%
    tidyr::separate(sample_id, c("date", "mouse", "tissue"), remove = FALSE) %>%
    mutate(short_label = sprintf("mouse %s, %s°C", mouse, condition_temp)) %>%
    mutate(long_label = paste(sample_id, condition_temp)) %>%
    mutate(colour = brewer.pal(length(unique(condition_temp)), "Set1")[condition_temp]) %>%
    arrange(condition_temp)

# create directories
dir.create(file.path(hub, "mm10"),
    showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(bw),
    showWarnings = FALSE)

# create genomes.txt
cat("genome mm10\ntrackDb mm10/trackDb.txt",
    file = file.path(hub, "genomes.txt"))

# create hub.txt
cat("hub nanoporeibat",
    "shortLabel illumina reads of murine iBAT",
    "longLabel RNA from BAT of cold treated and control mice.",
    "genomesFile genomes.txt",
    "email  christophak@bmb.sdu.dk",
    "descriptionUrl",
    sep = '\n',
    file = file.path(hub, "hub.txt"))

# create trackDb.txt and copy files

# remove if exists
suppressWarnings(file.remove(file.path(hub, "mm10", "trackDb.txt")))

# illumina
for (i in 1:nrow(sample_info)) {

    row <- sample_info[i,]

    for (direction in c("fw", "rv")) {

        file <- sprintf("%s_illumina_%s.bw", row$sample_id, direction)

        track <- c(
            sprintf("track %s_illumina_%s", row$sample_id, direction),
            sprintf("bigDataUrl %s", paste0(snakemake@params[["url"]], "/bw/", file)),
            sprintf("shortLabel %s, illumina %s", row$short_label, direction),
            sprintf("longLabel %s, illumina %s", row$long_label, direction),
            "type bigWig",
            "visibility full",
            "autoScale on",
            sprintf("color %s", col2rgb(row$colour) %>%
                as.vector() %>%
                paste(collapse = ',')),
            "")

        cat(track,
            file = file.path(hub, "mm10", "trackDb.txt"),
            sep = '\n',
            append = TRUE)

        file.copy(file.path("bw/illumina", file), file.path(bw, file))
    }
}

# ont
for (dataset in c("cdna", "teloprime")) {

    for (i in 1:nrow(sample_info)) {

        row <- sample_info[i,]

        file <- sprintf("%s_%s.bw", row$sample_id, dataset)

        track <- c(
            sprintf("track %s_%s", row$sample_id, dataset),
            sprintf("bigDataUrl %s", paste0(snakemake@params[["url"]], "/bw/", file)),
            sprintf("shortLabel %s, %s", row$short_label, dataset),
            sprintf("longLabel %s, %s", row$long_label, dataset),
            "type bigWig",
            "visibility full",
            "autoScale on",
            sprintf("color %s", col2rgb(row$colour) %>%
                as.vector() %>%
                paste(collapse = ',')),
            "")

        cat(track,
            file = file.path(hub, "mm10", "trackDb.txt"),
            sep = '\n',
            append = TRUE)

        file.copy(file.path("bw", dataset, file), file.path(bw, file))
    }
}

# flair
for (file in c("nanoporeibat_hub/bigGenePred/flair_teloprime.bb",
               "nanoporeibat_hub/bigGenePred/flair_cdna.bb")) {

    name <- basename(file)

    track <- c(
        sprintf("track %s", c(tools::file_path_sans_ext(name))),
        sprintf("bigDataUrl %s",
            file.path(snakemake@params[["url"]], "bigGenePred", name)),
        sprintf("shortLabel %s", tools::file_path_sans_ext(name)),
        sprintf("longLabel %s", tools::file_path_sans_ext(name)),
        "type bigGenePred",
        "visibility full",
        "")

        cat(track,
            file = file.path(hub, "mm10", "trackDb.txt"),
            sep = '\n',
            append = TRUE)
}
