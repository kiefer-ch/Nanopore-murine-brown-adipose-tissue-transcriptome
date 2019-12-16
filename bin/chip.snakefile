rule annotate_peaks:
    input:
        txdb = "annotation/annotation_txdb.sqlite",
        h3k4_warm = "data/chip/k4me3/NW_broad.bed",
        h3k4_cold = "data/chip/k4me3/NC_broad.bed"
    output:
        TSS = "res/chip/tss.csv.gz",
        full_range = "res/chip/full_range.csv.gz"
    script:
        "chip.R"
