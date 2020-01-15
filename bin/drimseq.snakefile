# illumina reads
rule drimseq_dmdsFromScaledTpm:
    input:
        txdb = "annotation/annotation_txdb.sqlite",
        tpm = "res/dexseq/illumina/dexseq_scaledTPM.rds",
        sample_info = "sample_info/sampleInfo.csv"
    output:
        "res/drimseq/illumina/illumina_dmds.rds"
    script:
        "drimseq_dmdsFromScaledTpm.R"


# ont
rule drimseq_dmdsFromCountMatrix_flair:
    input:
        txdb = "flair/{dataset}/flair.collapse.isoforms_txdb.sqlite",
        counts = "flair/{dataset}/flair_teloprime_counts_matrix.tsv",
        sample_info = "sample_info/sampleInfo.csv"
    output:
        "res/drimseq/{dataset}_flair/{dataset}_flair_dmds.rds"
    wildcard_constraints:
        dataset = "teloprime"
    script:
        "drimseq_dmdsFromCountMatrix_flair.R"

rule drimseq_dmdsFromCountMatrix:
    input:
        txdb = "annotation/annotation_txdb.sqlite",
        counts =
        sample_info = "sample_info/sampleInfo.csv"


# common
rule drimseq_stageR:
    input:
        dmds = "{file}_dmds.rds",
        biomaRt_gene = "annotation/biomaRt_gene.rds",
        biomaRt_tx = "annotation/biomaRt_tx.rds"
    params:
        out_folder = "{file}"
    output:
        "{file}_drimSeqStageR.html"
    threads:
        4
    script:
        "drimseq_stagerAnalysis.Rmd"


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


rule drimseq_browserPlots_flair:
    params:
        out_folder = "res/drimseq/teloprime_flair/browser_plots"
    input:
        genome = "annotation/genome.fa",
        bw_warm = "bw/teloprime/barcode01.bw",
        bw_cold = "bw/teloprime/barcode02.bw",
        txdb = "flair/teloprime/flair.collapse.isoforms_txdb.sqlite",
        stageR_results = "res/drimseq/teloprime_flair/teloprime_flair_drimSeqStageR.csv"
    script:
        "drimseq_browserPlots_flair.R"
