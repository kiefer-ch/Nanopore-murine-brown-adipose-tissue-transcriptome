#!/usr/bin/Rscript --no-restore --no-environ --no-save

# set libpaths to packrat local library
source("packrat/init.R")

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

save.image("browser_plots.RData")

txdb <- loadDb(snakemake@input[["txdb"]])
biomart <- read_rds(snakemake@input[["biomaRt_tx"]])

plot.gene <- function(gene = "ENSMUSG00000020654.15", extend_right = 50) {
    df <- AnnotationDbi::select(txdb, keys = gene,
            keytype = "GENEID",
            columns = c("TXCHROM", "TXSTART", "TXEND", "TXSTRAND")) %>%
        as_tibble() %>%
        group_by(GENEID) %>%
        summarise(TXCHROM = unique(TXCHROM),
            TXSTART = min(TXSTART),
            TXEND = max(TXEND),
            TXSTRAND = unique(TXSTRAND))
    strand <- ifelse(df$TXSTRAND == '+', "fw", "rv")

#    itrack <- IdeogramTrack(genome = "mm10", chromosome = df$TXCHROM)
    gtrack <- GenomeAxisTrack(genome = "mm10", chromosome = df$TXCHROM)
    atrack <- GeneRegionTrack(txdb, chromosome = df$TXCHROM,
        from = df$TXSTART, to = df$TXEND,
        name = "GENCODE M23")
    symbol(atrack) <- symbol(atrack) %>%
        tibble::enframe() %>%
        dplyr::select(ensembl_transcript_id_version = "value") %>%
        left_join(biomart, by = "ensembl_transcript_id_version") %>%
        pull(mgi_symbol)
    dtrack1 <- DataTrack(range = file.path(snakemake@params[["bw_folder"]], sprintf("5034_S33_%s.bw", strand)),
        type = "hist",
        genome = "mm10",
        name = "22°C",
        chromosome = df$TXCHROM,
        start = df$TXSTART,
        end = df$TXEND + extend_right)
    dtrack2 <- DataTrack(range = file.path(snakemake@params[["bw_folder"]], sprintf("5035_S34_%s.bw", strand)),
        type = "hist",
        genome = "mm10",
        name = "4°C",
        chromosome = df$TXCHROM,
        start = df$TXSTART,
        end = df$TXEND + extend_right)
    # dtrack3 <- DataTrack(range = file.path(snakemake@params[["bw_folder"]], sprintf("5039_S38_%s.bw", strand)),
    #     type = "hist",
    #     genome = "mm10",
    #     name = "warm2",
    #     chromosome = df$TXCHROM,
    #     start = df$TXSTART,
    #     end = df$TXEND)
    # dtrack4 <- DataTrack(range = file.path(snakemake@params[["bw_folder"]], sprintf("5041_S40_%s.bw", strand)),
    #     type = "hist",
    #     genome = "mm10",
    #     name = "cold2",
    #     chromosome = df$TXCHROM,
    #     start = df$TXSTART,
    #     end = df$TXEND)
    # dtrack5 <- DataTrack(range = file.path(snakemake@params[["bw_folder"]], sprintf("5043_S42_%s.bw", strand)),
    #     type = "hist",
    #     genome = "mm10",
    #     chromosome = df$TXCHROM,
    #     start = df$TXSTART,
    #     end = df$TXEND)
    # dtrack6 <- DataTrack(range = file.path(snakemake@params[["bw_folder"]], sprintf("5044_S43_%s.bw", strand)),
    #     type = "hist",
    #     genome = "mm10",
    #     chromosome = df$TXCHROM,
    #     start = df$TXSTART,
    #     end = df$TXEND)

    plotTracks(list(
#            itrack,
            gtrack,
            dtrack1, #dtrack3, #dtrack5,
            dtrack2, #dtrack4, #dtrack6,
            atrack),
#        showId = TRUE, geneSymbol = TRUE,
        transcriptAnnotation = "symbol",
        from = df$TXSTART,
        to = df$TXEND + extend_right,
        extend.left = 50,
        extend_right = 50,
        background.title = "darkblue",
#        collapseTranscripts = "meta",
        fontsize = 15, cex.axis = .4)
}

plot.bam <- function(gene = "ENSMUSG00000020654.15", extend_right = 50,
        type = "both") {
    df <- AnnotationDbi::select(txdb, keys = gene,
            keytype = "GENEID",
            columns = c("TXCHROM", "TXSTART", "TXEND", "TXSTRAND")) %>%
        as_tibble() %>%
        group_by(GENEID) %>%
        summarise(TXCHROM = unique(TXCHROM),
            TXSTART = min(TXSTART),
            TXEND = max(TXEND),
            TXSTRAND = unique(TXSTRAND))
    strand <- ifelse(df$TXSTRAND == '+', "fw", "rv")

    sTrack <- SequenceTrack("annotation/genome.fa",
        genome = "mm10", chromosome = df$TXCHROM)

    gtrack <- GenomeAxisTrack(genome = "mm10", chromosome = df$TXCHROM)
    atrack <- GeneRegionTrack(txdb, chromosome = df$TXCHROM,
        from = df$TXSTART, to = df$TXEND, strand = df$TXSTRAND,
        name = "GENCODE M23")
    symbol(atrack) <- symbol(atrack) %>%
        tibble::enframe() %>%
        dplyr::select(ensembl_transcript_id_version = "value") %>%
        left_join(biomart, by = "ensembl_transcript_id_version") %>%
        pull(mgi_symbol)

    if (type %in% c("illumina", "both")) {
        illumina <- AlignmentsTrack(snakemake@input[["illumina_bam"]],
            isPaired = TRUE,
            genome = "mm10",
            name = "illumina",
            chromosome = df$TXCHROM,
            start = df$TXSTART,
            end = df$TXEND,
            referenceSequence = sTrack)
    }

    if (type %in% c("teloprime", "both")) {
        teloprime <- AlignmentsTrack(snakemake@input[["teloprime_bam"]],
            isPaired = FALSE,
            genome = "mm10",
            name = "telorprime",
            chromosome = df$TXCHROM,
            start = df$TXSTART,
            end = df$TXEND,
            referenceSequence = sTrack)
    }

    options(ucscChromosomeNames=FALSE)

    plotTracks(list(
            gtrack,
    #        illumina,
            teloprime,
            atrack),
        transcriptAnnotation = "symbol",
        from = df$TXSTART,
        to = df$TXEND,
        extend.left = 50,
        extend.right = extend_right,
        background.title = "darkblue",
        fontsize = 15, cex.axis = .4,
        min.height = 0, coverageHeight = 0.08, minCoverageHeight = 0)
}


pdf("adcy3.pdf", width = 12, height = 12 / 1.618)
    plot.gene()
dev.off()

pdf("adcy3_bam.pdf", width = 12, height = 12 / 1.618)
    plot.bam(type = "teloprime")
dev.off()

pdf("ctcflos.pdf", width = 12, height = 12 / 1.618)
    plot.gene(gene = "ENSMUSG00000087382.7",  extend_right = 6000)
dev.off()

pdf("ctcflos_bam.pdf", width = 12, height = 12 / 1.618)
    plot.bam(gene = "ENSMUSG00000087382.7", extend_right = 6000)
dev.off()

pdf("Gm15551.pdf", width = 12, height = 12 / 1.618)
    plot.gene(gene = "ENSMUSG00000086679.2")
dev.off()

pdf("Gm15551.pdf", width = 12, height = 12 / 1.618)
    plot.bam(gene = "ENSMUSG00000086679.2")
dev.off()
