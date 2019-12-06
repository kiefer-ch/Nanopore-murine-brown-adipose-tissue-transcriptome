# DTU
rule drimseq_ont_dtu:
    threads: 4
    params:
        out_folder = "res/drimseq/dtu_ont"
    input:
        annotation = "annotation/annotation.gtf",
        tpm = "data/scaledTPM_ont.rds",
        sample_info = "sample_info/sampleInfo.csv"
    output:
        "res/drimseq/dtu_ont/ont_dtu.html"
    script:
        "drimseq_ont.Rmd"

rule drimseq_all_dtu:
    threads: 4
    params:
        out_folder = "res/drimseq/dtu_all"
    input:
        annotation = "annotation/annotation.gtf",
        tpm = "data/scaledTPM_all.rds",
        sample_info = "sample_info/sampleInfo.csv"
    output:
        "res/drimseq/dtu_all/all_dtu.html"
    script:
        "drimseq_all.Rmd"

rule drimseq_dmdsFromCountMatrix:
    input:
        txdb = "flair/{dataset}/flair.collapse.isoforms_txdb.sqlite",
        counts = "flair/{dataset}/flair_teloprime_counts_matrix.tsv",
        sample_info = "sample_info/sampleInfo.csv"
    output:
        "res/drimseq/{dataset}_flair/{dataset}_flair_dmds.rds"
    wildcard_constraints:
        dataset = "teloprime"
    script:
        "drimseq_dmdsFromCountMatrix.R"

rule drimseq_browserPlots:
    params:
        out_folder = "res/drimseq/teloprime/browser_plots"
    input:
        genome = "annotation/genome.fa",
        bw_warm = "bw/teloprime/barcode01.bw",
        bw_cold = "bw/teloprime/barcode02.bw",
        txdb = "annotation/annotation_txdb.sqlite",
        biomaRt_tx = "annotation/biomaRt_tx.rds",
        stageR_results = "res/wien/teloprime_old/DRIMSeq_stageR/stageR/stageR_final_output_padj_GeneSymbols.tsv"
    script:
        "drimseq_browserPlots.R"
