# DTU
rule drimseq_ont_dtu:
    threads: 4
    params:
        out_folder = "res/dtu_ont"
    input:
        annotation = "annotation/annotation.gtf",
        tpm = "data/scaledTPM_ont.rds",
        sample_info = "sample_info/sampleInfo.csv"
    output:
        "res/dtu/dtu_ont/ont_dtu.html"
    script:
        "drimseq_ont.Rmd"

rule drimseq_all_dtu:
    threads: 4
    params:
        out_folder = "res/dtu_all"
    input:
        annotation = "annotation/annotation.gtf",
        tpm = "data/scaledTPM_all.rds",
        sample_info = "sample_info/sampleInfo.csv"
    output:
        "res/dtu/dtu_all/all_dtu.html"
    script:
        "drimseq_all.Rmd"
