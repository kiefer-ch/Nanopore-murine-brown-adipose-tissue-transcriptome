rule browserTrackAdcy3:
    input:
        "data/bam/illumina/5034_S33_Aligned.sortedByCoord.out.bam.bai",
        "data/bam/illumina/5035_S34_Aligned.sortedByCoord.out.bam.bai",
        "data/bam/cdna/merged/cdna_barcode07_genome.bam.bai",
        "data/bam/cdna/merged/cdna_barcode08_genome.bam.bai",
        txdb = "data/annotation/annotation_txdb.sqlite",
        txdb_flair = "data/reannotation/flair/annotation/cdna_flair.isoforms_txdb.sqlite",
        txdb_stringtie = "data/reannotation/stringtie/teloprime_stringtie_noUnknownStrand_txdb.sqlite",
        tmap = "data/reannotation/stringtie/gffcmp.teloprime_stringtie.gtf.tmap",
        txdb_stringtie_illumina = "data/reannotation/stringtie/illumina_stringtie_noUnknownStrand_txdb.sqlite",
        tmap_illumina = "data/reannotation/stringtie/gffcmp.illumina_stringtie.gtf.tmap",
        illumina_warm = "data/bam/illumina/5034_S33_Aligned.sortedByCoord.out.bam",
        illumina_cold = "data/bam/illumina/5035_S34_Aligned.sortedByCoord.out.bam",
        ont_warm = "data/bam/cdna/merged/cdna_barcode07_genome.bam",
        ont_cold = "data/bam/cdna/merged/cdna_barcode08_genome.bam",
        ncd_warm_h3k4me3 = "data/additionalBrowserTracks/200901_rute_chip_adipocytes/2_warm_ncd_H3K4me3.bw",
        ncd_cold_h3k4me3 =  "data/additionalBrowserTracks/200901_rute_chip_adipocytes/7_warm_hfd_H3K4me3.bw",
        genome = "data/annotation/genome.fa.fai"
    output:
        "res/browserTracks/Adcy3.pdf"
    params:
        gene_id = "ENSMUSG00000020654.15",
        max_cov_illumina = 5,
        lwd_sashimi_max = 10,
        extend_plot = 0,
        temp = "warm",
        reverse_strand = False,
        dimension = [6.85, 6]
    conda:
        "../envs/r_4.1.2.yaml"
    script:
        "browserTracks.R"


rule browserTrackCars2:
    input:
        "data/bam/illumina/5034_S33_Aligned.sortedByCoord.out.bam.bai",
        "data/bam/illumina/5035_S34_Aligned.sortedByCoord.out.bam.bai",
        "data/bam/cdna/merged/cdna_barcode07_genome.bam.bai",
        "data/bam/cdna/merged/cdna_barcode08_genome.bam.bai",
        txdb = "data/annotation/annotation_txdb.sqlite",
        txdb_flair = "data/reannotation/flair/annotation/cdna_flair.isoforms_txdb.sqlite",
        txdb_stringtie = "data/reannotation/stringtie/teloprime_stringtie_noUnknownStrand_txdb.sqlite",
        tmap = "data/reannotation/stringtie/gffcmp.teloprime_stringtie.gtf.tmap",
        txdb_stringtie_illumina = "data/reannotation/stringtie/illumina_stringtie_noUnknownStrand_txdb.sqlite",
        tmap_illumina = "data/reannotation/stringtie/gffcmp.illumina_stringtie.gtf.tmap",
        illumina_warm = "data/bam/illumina/5034_S33_Aligned.sortedByCoord.out.bam",
        illumina_cold = "data/bam/illumina/5035_S34_Aligned.sortedByCoord.out.bam",
        ont_warm = "data/bam/cdna/merged/cdna_barcode07_genome.bam",
        ont_cold = "data/bam/cdna/merged/cdna_barcode08_genome.bam",
        ncd_warm_h3k4me3 = "data/additionalBrowserTracks/200901_rute_chip_adipocytes/2_warm_ncd_H3K4me3.bw",
        ncd_cold_h3k4me3 =  "data/additionalBrowserTracks/200901_rute_chip_adipocytes/7_warm_hfd_H3K4me3.bw",
        genome = "data/annotation/genome.fa.fai"
    output:
        "res/browserTracks/Cars2.pdf"
    params:
        gene_id = "ENSMUSG00000056228.10",
        max_cov_illumina = 6,
        lwd_sashimi_max = 10,
        temp = "both",
        extend_plot = 400,
        reverse_strand = True,
        dimension = [6.85, 6]
    conda:
        "../envs/r_4.1.2.yaml"
    script:
        "browserTracks.R"


rule makeBrowserTracks:
    input:
        "res/browserTracks/Adcy3.pdf",
        "res/browserTracks/Cars2.pdf"
