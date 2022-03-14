rule qpcr_temperature_effect:
    input:
        cq = "data/qpcr/190514/190514_bl6_coldIBat.txt",
        sample_info = config["SAMPLE_INFO"]
    output:
        "res/qpcr/qpcr_temperature_effect.html"
    script:
        "qpcr_temperature_effect.Rmd"


rule dtu_qpcr:
    input:
        cq1 = "data/qpcr/200514/200514_bl6_dtu_1.txt",
        cq2 = "data/qpcr/200514/200514_bl6_dtu_2.txt",
        cq3 = "data/qpcr/200611/200611_bl6_dtu_3.txt",
        sample_info = config["SAMPLE_INFO"]
    output:
        "data/qpcr/dtu_qpcr.csv",
        "data/qpcr/dtu_qpcr_cq.rds"
    script:
        "qpcr_dtu.R"


rule qpcr_dtu_report:
    input:
        cq = "data/qpcr/dtu_qpcr_cq.rds",
        txdb = "data/annotation/annotation_txdb.sqlite",
        txdb_flair = "data/reannotation/flair/annotation/cdna_flair.isoforms_txdb.sqlite",
        txdb_stringtie = "data/reannotation/stringtie/teloprime_stringtie_noUnknownStrand_txdb.sqlite",
        tmap = "data/reannotation/stringtie/gffcmp.teloprime_stringtie.gtf.tmap",
        txdb_stringtie_illumina = "data/reannotation/stringtie/illumina_stringtie_noUnknownStrand_txdb.sqlite",
        tmap_illumina = "data/reannotation/stringtie/gffcmp.illumina_stringtie.gtf.tmap",
#        illumina_warm = "data/bam/illumina/5034_S33_Aligned.sortedByCoord.out.bam",
#        illumina_cold = "data/bam/illumina/5035_S34_Aligned.sortedByCoord.out.bam",
#        cdna_warm = "data/bam/cdna/merged/cdna_barcode07_genome.bam",
#        cdna_cold = "data/bam/cdna/merged/cdna_barcode08_genome.bam",
#        teloprime_warm = "data/bam/cdna/merged/cdna_barcode01_genome.bam",
#        teloprime_cold = "data/bam/cdna/merged/cdna_barcode02_genome.bam",
        genome = "data/annotation/genome.fa.fai",
        primer_stats = "data/qpcr/200331_primer.csv"
    output:
        "res/qpcr/qpcr_dtu_validation.html"
    script:
        "qpcr_dtu_report.Rmd"