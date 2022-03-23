
source(".Rprofile")
suppressPackageStartupMessages({
    library("GenomicFeatures")
    library("Gviz")
    library("dplyr")
    library("readr")
    library("logger")
    library("AnnotationDbi")
    source("R/subsetTxdb.R")
})


log_info("Reading annotation files..")
txdb_gencode <- loadDb(snakemake@input[["txdb"]])
txdb_gencode_dump <- as.list(txdb_gencode)

txdb_flair <- loadDb(snakemake@input[["txdb_flair"]])
txdb_flair_dump <- as.list(txdb_flair)

txdb_stringtie <- loadDb(snakemake@input[["txdb_stringtie"]])
txdb_stringtie_dump <- as.list(txdb_stringtie)

stringtie_refGene <- read_tsv(snakemake@input$tmap,
    col_types = "cccccidddici",
    na = "-") %>%
    dplyr::select(ref_gene_id, qry_gene_id) %>%
    tidyr::drop_na() %>%
    distinct()

txdb_stringtie_illumina <- loadDb(snakemake@input[["txdb_stringtie_illumina"]])
txdb_stringtie_illumina_dump <- as.list(txdb_stringtie_illumina)

stringtie_illumina_refGene <- read_tsv(snakemake@input$tmap_illumina,
                              col_types = "cccccidddici",
                              na = "-") %>%
    dplyr::select(ref_gene_id, qry_gene_id) %>%
    tidyr::drop_na() %>%
    distinct()


faidx <- read_tsv(snakemake@input$genome,
        col_types = "cidii",
        col_names = c("chrom", "length", "offset", "LINEBASES", "LINEWIDTH")) %>%
    dplyr::select(chrom, length) %>%
    mutate(is_circular = FALSE) %>%
    as.data.frame()



# function definitions
plot.transcripts <- function(gene_id, max_cov_illumina = 5,
    extend = 2500, lwd_sashimi_max = 3, temp = "warm") {

    gene_info <- suppressMessages(AnnotationDbi::select(txdb_gencode, keys = gene_id,
            keytype = "GENEID",
            columns = c("TXCHROM", "TXSTART", "TXEND"))) %>%
        as_tibble() %>%
        group_by(GENEID) %>%
        summarise(TXCHROM = unique(TXCHROM),
            TXSTART = min(TXSTART),
            TXEND = max(TXEND))


    # gencode
    txdb_gencode_subset <- subset.txdb(txdb_gencode, txdb_gencode_dump, gene_id)

    gencodetrack <- GeneRegionTrack(txdb_gencode_subset,
        chromosome = gene_info$TXCHROM,
        min.height = 3,
        from = gene_info$TXSTART, to = gene_info$TXEND,
        name = "M22")


    # flair
    txdb_flair_subset <- subset.txdb(txdb_flair, txdb_flair_dump, gene_id)

    flairtrack <- GeneRegionTrack(txdb_flair_subset,
        chromosome = gene_info$TXCHROM,
        min.height = 3,
        from = gene_info$TXSTART, to = gene_info$TXEND,
        name = "flair cDNA")


    # stringtie
    txdb_stringtie_subset <- subset.txdb(txdb_stringtie, txdb_stringtie_dump,
        stringtie_refGene %>%  filter(ref_gene_id == gene_id) %>%  pull(qry_gene_id))

    stringtietrack <- GeneRegionTrack(txdb_stringtie_subset,
        chromosome = gene_info$TXCHROM,
        min.height = 3,
        from = gene_info$TXSTART, to = gene_info$TXEND,
        name = "stringtie TeloPrime")

    # stringtie illumina
    txdb_stringtie_illumina_subset <- subset.txdb(txdb_stringtie_illumina, txdb_stringtie_illumina_dump,
        stringtie_illumina_refGene %>%  filter(ref_gene_id == gene_id) %>%  pull(qry_gene_id))

    stringtietrack2 <- GeneRegionTrack(txdb_stringtie_illumina_subset,
                                      chromosome = gene_info$TXCHROM,
                                      min.height = 3,
                                      from = gene_info$TXSTART, to = gene_info$TXEND,
                                      name = "stringtie Illumina")



    # list of unique exons for sashimi plot
    unique_exons <- c(txdb_flair_subset, txdb_stringtie_subset, txdb_gencode_subset) %>%
        purrr::map(intronsByTranscript) %>% # extra intron information
        purrr::reduce(c) %>% # collapse the lists into one single list
        unlist() %>% # concatenate all exons of all transcripts into one granges object
        unique() # remove duplicate exons


    # scale
    gtrack <- GenomeAxisTrack(from = min(gene_info$TXSTART),
        to = max(gene_info$TXEND),
        scale = 0.25)


    # bam tracks
    if(temp == "warm") {

        ont <- AlignmentsTrack(snakemake@input[["ont_warm"]],
            isPaired = FALSE,
            genome = "mm10",
            name = "cDNA 22°C",
            type = c("pileup"),
            min.height = 3, # minimum height of a read in px, controls how many read are plotted
            chromosome = gene_info$TXCHROM,
            start = gene_info$TXSTART,
            end = gene_info$TXEND)

    } else if (temp == "cold") {

        ont <- AlignmentsTrack(snakemake@input[["ont_cold"]],
            isPaired = FALSE,
            genome = "mm10",
            name = "cDNA 4°C",
            type = c("pileup"),
            min.height = 3, # minimum height of a read in px, controls how many read are plotted
            chromosome = gene_info$TXCHROM,
            start = gene_info$TXSTART,
            end = gene_info$TXEND)

    }

    # illumina

    illumina_w <- AlignmentsTrack(snakemake@input[["illumina_warm"]],
                                isPaired = TRUE,
                                genome = "mm10",
                                name = "Illumina 22°C",
                                type = c("coverage", "sashimi"),
                                ylim = c(0, max_cov_illumina),
                                transformation = function(x) log1p(x),
                                sashimiTransformation = function(x) x,
                                lwd.sashimiMax = lwd_sashimi_max,
                                sashimiHeight = .33,
                                chromosome = gene_info$TXCHROM,
                                start = gene_info$TXSTART,
                                end = gene_info$TXEND)

    illumina_c <- AlignmentsTrack(snakemake@input[["illumina_cold"]],
                                isPaired = TRUE,
                                genome = "mm10",
                                name = "Illumina 4°C",
                                type = c("coverage", "sashimi"),
                                ylim = c(0, max_cov_illumina),
                                transformation = function(x) log1p(x),
                                sashimiTransformation = function(x) x,
                                lwd.sashimiMax = lwd_sashimi_max,
                                sashimiHeight = .33,
                                chromosome = gene_info$TXCHROM,
                                start = gene_info$TXSTART,
                                end = gene_info$TXEND)


    # chip
    h3k4me3 <- DataTrack(snakemake@input[["ncd_cold_h3k4me3"]],
        genome = "mm10",
        name = "H3K4me3",
        chromosome = unique(gene_info$TXCHROM),
        start = min(gene_info$TXSTART) - extend,
        end = max(gene_info$TXEND) + extend,
        type = "histogram")


    # Plotting
    options(ucscChromosomeNames = FALSE)

    tracks <- list(gtrack,
        illumina_w, illumina_c,
        ont,
        gencodetrack, stringtietrack2, stringtietrack, flairtrack,
        h3k4me3)


    if(length(extend) == 1) {

        extend_l <- extend
        extend_r <- extend

    } else if (length(extend) == 2) {

        extend_l <- extend[1]
        extend_r <- extend[2]

    }

    plotTracks(
        tracks,
        from = min(gene_info$TXSTART) - extend_l,
        to = max(gene_info$TXEND) + extend_r,
        extend.left = 250,
        extend.right = 250,
        background.title = "transparent",
        col.title = "grey10",
        cex = .75,
        cex.axis = .6,
        sashimiFilter = unique_exons,
        lwd = .25,
        sizes = c(.1, rep(.33, 2), 1, rep(.25, 4), .25))

}

log_info(sprintf("Plotting %s...", snakemake@params$gene_id))

pdf(snakemake@output[[1]],
    width = snakemake@params$dimsenion[1],
    height = snakemake@params$dimension[2])

    plot.transcripts(
        gene_id = snakemake@params$gene_id,
        max_cov_illumina = snakemake@params$max_cov_illumina,
        lwd_sashimi_max = snakemake@params$lwd_sashimi_max,
        extend = snakemake@params$extend_plot,
        temp = snakemake@params$temp)

dev.off()
