rule browserTracksAdcy3Cars2:
    input:
        txdb = "annotation/annotation_txdb.sqlite",
        txdb_flair = "flair/cdna/flair.collapse.cdna.isoforms_txdb.sqlite",
        illumina_warm = "bam/illumina/5034_S33_Aligned.sortedByCoord.out.bam",
        illumina_cold = "bam/illumina/5035_S34_Aligned.sortedByCoord.out.bam",
        cdna_warm = "bam/cdna/barcode07_genome.bam",
        cdna_cold = "bam/cdna/barcode08_genome.bam",
        ncd_warm_h3k4me3 = "data/additionalBrowserTracks/200901_rute_chip_adipocytes/2_warm_ncd_H3K4me3.bw",
        ncd_cold_h3k4me3 =  "data/additionalBrowserTracks/200901_rute_chip_adipocytes/7_warm_hfd_H3K4me3.bw",
        genome = "annotation/genome.fa.fai"
    output:
        "res/browserTracks/nanopore_iBAT_browserTracks_Adcy3Cars2.html"
    script:
        "browserTracks_Adcy3Cars2.Rmd"
