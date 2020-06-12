rule qpcr_temperature_effect:
    input:
        cq = "data/qpcr/190514/190514_bl6_coldIBat.txt",
        sample_info = "sample_info/sampleInfo.csv"
    output:
        "res/qpcr/qpcr_temperature_effect.html"
    script:
        "qpcr_temperature_effect.Rmd"

rule qpcr_dtu_validation:
    input:
        cq1 = "data/qpcr/200514/200514_bl6_dtu_1.txt",
        cq2 = "data/qpcr/200514/200514_bl6_dtu_2.txt",
        cq3 = "data/qpcr/200611/200611_bl6_dtu_3.txt",
        sample_info = "sample_info/sampleInfo.csv",
        signif = "res/comparisons/comparisons_dtu_all.rds",
        txdb = "annotation/annotation_txdb.sqlite",
        txdb_flair = "flair/cdna/flair.collapse.isoforms_txdb.sqlite",
        illumina_warm = "bam/illumina/5034_S33_Aligned.sortedByCoord.out.bam",
        illumina_cold = "bam/illumina/5035_S34_Aligned.sortedByCoord.out.bam",
        teloprime_warm = "bam/teloprime/barcode01_genome.bam",
        teloprime_cold = "bam/teloprime/barcode02_genome.bam",
        cdna_warm = "bam/cdna/barcode07_genome.bam",
        cdna_cold = "bam/cdna/barcode08_genome.bam",
        genome = "annotation/genome.fa",
        primer_stats = "res/christoph/primer_design/200331_primer.csv"
    output:
        "res/qpcr/qpcr_dtu_validation.html"
    script:
        "qpcr_dtu_validation.Rmd"
