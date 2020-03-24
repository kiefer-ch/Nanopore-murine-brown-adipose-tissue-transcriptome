#!/usr/bin/Rscript --no-restore --no-environ --no-save

# set libpaths to packrat local library
source(".Rprofile")

suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("Gviz"))
suppressPackageStartupMessages(library("Mus.musculus"))

################################################################################
#
# author: Christoph Kiefer
# email: christophak@bmb.sdu.dk
#
################################################################################

stageR <- read_csv(snakemake@input[["stageR_results"]])
txdb <- loadDb(snakemake@input[["txdb"]])

plot.bw <- function(gene_id = "ENSMUSG00000027327.16") {
    # stageR significant transcripts
    sig_tx <- stageR %>%
        filter(ensembl_gene_id_version == gene_id) %>%
        filter(transcript < .05) %>%
        pull(ensembl_transcript_id_version)

    # df with chromosome, tx_start, tx_end and strand
    df <- AnnotationDbi::select(txdb, keys = gene_id,
            keytype = "GENEID",
            columns = c("TXCHROM", "TXSTART", "TXEND", "TXSTRAND")) %>%
        as_tibble() %>%
        group_by(GENEID) %>%
        summarise(TXCHROM = unique(TXCHROM),
            TXSTART = min(TXSTART),
            TXEND = max(TXEND),
            TXSTRAND = unique(TXSTRAND))
    strand <- ifelse(df$TXSTRAND == '+', "fw", "rv")

    # genome track
    sTrack <- SequenceTrack(snakemake@input[["genome"]],
        genome = "mm10", chromosome = df$TXCHROM)

    # annotation track
    gtrack <- GenomeAxisTrack(genome = "mm10", chromosome = df$TXCHROM)
    atrack <- GeneRegionTrack(txdb, chromosome = df$TXCHROM,
        from = df$TXSTART, to = df$TXEND, strand = df$TXSTRAND,
        name = "GENCODE M23")

    # set colour of significant transcripts
    feature(atrack) <- if_else(transcript(atrack) %in% sig_tx, "bad", "good")
    interestcolor <- list("bad" = "red", "good" = "#FFD58A")
    displayPars(atrack) <- interestcolor

    # bam tracks
    bw_warm <- DataTrack(range = snakemake@input[["bw_warm"]],
        type = "hist",
        genome = "mm10",
        name = "22°C",
        chromosome = df$TXCHROM,
        start = df$TXSTART,
        end = df$TXEND)
    bw_cold <- DataTrack(range = snakemake@input[["bw_cold"]],
        type = "hist",
        genome = "mm10",
        name = "4°C",
        chromosome = df$TXCHROM,
        start = df$TXSTART,
        end = df$TXEND)

    options(ucscChromosomeNames=FALSE)

    tracklist <- list(gtrack, atrack, bw_warm, bw_cold)

    plotTracks(tracklist,
        just.group="above",
        transcriptAnnotation = "transcript",
        from = df$TXSTART,
        to = df$TXEND,
        extend.left = 250,
        extend.right = 250,
        background.title = "darkblue",
        fontsize = 15, cex.axis = .4,
        min.height = 0, coverageHeight = 0.08, minCoverageHeight = 0)
}

genes <- stageR %>%
    filter(gene < .05) %>%
    pull(ensembl_gene_id_version) %>%
    unique()

dir.create(snakemake@params[["out_folder"]],
    showWarnings = FALSE, recursive = TRUE)

for (gene in genes) {
    pdf(file.path(snakemake@params[["out_folder"]], paste0(gene, ".pdf")),
            height = 7, width = 7 * 1.618)
        try(plot.bw(gene_id = gene))
    dev.off()
}
