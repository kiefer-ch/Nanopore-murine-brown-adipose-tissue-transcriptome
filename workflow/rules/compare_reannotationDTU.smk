# reannotation
rule compare_reannotation:
    input:
        gffcompare =
            [expand("data/reannotation/flair/annotation/gffcmp.{dataset}_flair.isoforms.gtf.tmap",
                dataset=["cdna", "teloprime", "rna"]),
            expand("data/reannotation/stringtie/gffcmp.{dataset}_stringtie.gtf.tmap",
                dataset=["cdna", "teloprime", "rna", "illumina"])],
        gffcmp_tracking = "data/comparisons/reannotation/gffcompare/R/gffcmp.tracking",
        sqanti =
            [expand("data/comparisons/reannotation/squanti/stringtie/{dataset}/{dataset}_stringtie_noUnknownStrand_classification.txt",
                dataset=["cdna", "teloprime", "rna", "illumina"]),
            expand("data/comparisons/reannotation/squanti/flair/{dataset}/{dataset}_flair.isoforms_classification.txt",
                dataset=["cdna", "teloprime", "rna"])],
        gffcompare_stats = "data/comparisons/reannotation/gffcompare/RQ/gffcmp.stats",
        sal = expand("data/alternativeSpliceAnalysis/{dataset}_sal.rds",
            dataset = ["cdna_flair", "cdna_stringtie", "illumina_stringtie", "teloprime_flair",
                "teloprime_stringtie", "rna_flair", "rna_stringtie"])
    conda:
        "../envs/r_4.1.2.yaml"
    output:
        "res/comparisons/comparisons_reannotation.html"
    script:
        "../scripts/comparisons_reannotation.Rmd"


rule compare_dtu:
    input:
        dtu_res = expand("data/drimseq/{dataset}_{annotation}/{dataset}_{annotation}_dtu_res.csv.gz",
            dataset = ["illumina", "cdna"],
            annotation = ["ref", "illumina_stringtie", "cdna_flair", "teloprime_stringtie"]),
        sal = expand("data/drimseq/{dataset}_{annotation}/{dataset}_{annotation}_sal.rds",
            dataset = ["illumina", "cdna"],
            annotation = ["ref", "illumina_stringtie", "cdna_flair", "teloprime_stringtie"]),
        qpcr = "data/qpcr/dtu_qpcr.csv",
        annotation = "data/annotation/annotation_txdb.sqlite",
        telo_tmap = "data/reannotation/stringtie/gffcmp.teloprime_stringtie.gtf.tmap",
        illu_tmap = "data/reannotation/stringtie/gffcmp.illumina_stringtie.gtf.tmap"
    output:
        "res/comparisons/comparisons_dtu.html",
    params:
        dtu_cutoff = config["dtu_cutoff"], # alpha
        dIF_cutoff = config["dIF_cutoff"]  # min difference
    conda:
        "../envs/r_4.1.2.yaml"
    script:
        "../scripts/comparisons_dtu.Rmd"
