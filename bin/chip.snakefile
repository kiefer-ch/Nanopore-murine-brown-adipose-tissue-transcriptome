rule annotate_peaks_gencode:
    input:
        txdb = "annotation/annotation_txdb.sqlite",
        h3k4_warm = "data/chip/k4me3/NW_broad.bed",
        h3k4_cold = "data/chip/k4me3/NC_broad.bed"
    output:
        TSS = "res/chip/gencode_tss.csv.gz",
        full_range = "res/chip/gencode_full_range.csv.gz"
    script:
        "chip_annotatePeaks.R"

rule annotate_peaks_flair:
    input:
        txdb = "flair/{dataset}/flair.collapse.isoforms_txdb.sqlite",
        h3k4_warm = "data/chip/k4me3/NW_broad.bed",
        h3k4_cold = "data/chip/k4me3/NC_broad.bed"
    output:
        TSS = "res/chip/{dataset}_tss.csv.gz",
        full_range = "res/chip/{dataset}_full_range.csv.gz"
    wildcard_constraints:
        dataset = "teloprime"
    script:
        "chip_annotatePeaks.R"
