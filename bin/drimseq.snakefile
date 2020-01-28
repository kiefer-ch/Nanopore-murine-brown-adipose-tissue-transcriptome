# illumina reads
rule tximport_drimseq_illumina:
    input:
        salmon_out = expand("salmon/{sample}/quant.sf", sample=SAMPLES_ont),
        txdb = "annotation/annotation_txdb.sqlite",
        sample_info = "sample_info/sampleInfo.csv"
    output:
        "res/drimseq/illumina/drimseq_dtuScaledTPM.rds"
    script:
        "drimseq_txImport.R"


rule drimseq_dmdsFromScaledTpm:
    input:
        tpm = "res/drimseq/illumina/drimseq_dtuScaledTPM.rds",
        txdb = "annotation/annotation_txdb.sqlite",
        sample_info = "sample_info/sampleInfo.csv"
    output:
        "res/drimseq/illumina/illumina_dmds.rds"
    script:
        "drimseq_dmdsFromScaledTpm.R"


# ont
rule drimseq_dmdsFromCountMatrix_flair:
    input:
        txdb = "flair/{dataset}/flair.collapse.isoforms_txdb.sqlite",
        counts = "flair/{dataset}/flair_{dataset}_counts_matrix.tsv",
        sample_info = "sample_info/sampleInfo.csv"
    output:
        "res/drimseq/{dataset}_flair/{dataset}_flair_dmds.rds"
    wildcard_constraints:
        dataset = "teloprime|cdna"
    script:
        "drimseq_dmdsFromCountMatrix_flair.R"


rule drimseq_dmdsFromCountMatrix:
    input:
        txdb = "annotation/annotation_txdb.sqlite",
        counts = "res/deseq/{dataset}/txlevel/{dataset}_txlevel_cm_cts.csv.gz",
        sample_info = "sample_info/sampleInfo.csv"
    output:
        "res/drimseq/{dataset}/{dataset}_dmds.rds"
    wildcard_constraints:
        dataset = "teloprime|cdna"
    script:
        "drimseq_dmdsFromCountMatrix.R"


# common
rule drimseq_dtu:
    input:
        dmds = "{file}_dmds.rds",
        biomaRt_gene = "annotation/biomaRt_gene.rds"
    output:
        dmds = "{file}_dmds_dtu.rds",
        res = "{file}_drimSeqStageR.csv"
    threads:
        4
    script:
        "drimseq_dtu.R"


rule drimseq_report:
    input:
        dmds = "{file}_dmds_dtu.rds",
        res = "{file}_drimSeqStageR.csv",
        biomaRt_gene = "annotation/biomaRt_gene.rds",
        biomaRt_tx = "annotation/biomaRt_tx.rds"
    output:
        "{file}_drimSeqStageR.html"
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
        stageR_results = "res/drimseq/teloprime/teloprime_drimSeqStageR.csv"
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
