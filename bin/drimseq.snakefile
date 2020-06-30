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
        sample_info = "sample_info/sampleInfo.csv",
        dmds = "res/drimseq/{dataset}/{dataset}_dmds.rds",
        biomaRt_gene = "annotation/biomaRt_gene.rds"
    output:
        dmds = "res/drimseq/{dataset}/{dataset}_dmds_dtu.rds",
        res = "res/drimseq/{dataset}/{dataset}_drimSeqStageR.csv"
    threads:
        4
    script:
        "drimseq_dtu.R"


rule isoformswitchanalyser_importData:
    input:
        counts = "res/deseq/{dataset}/txlevel/{dataset}_txlevel_cm_cts.csv.gz",
        sample_info = "sample_info/sampleInfo.csv",
        gtf = "annotation/annotation.gtf",
        transcripts = "annotation/transcripts.fa",
        test = "res/drimseq/{dataset}/{dataset}_drimSeqStageR.csv"
    output:
        "res/drimseq/{dataset}/{dataset}_sal.rds"
    wildcard_constraints:
        dataset = "cdna|teloprime|illumina"
    script:
        "drimseq_isoformswitchanalyser_importData.R"


rule clean_fasta_ids:
    input:
        "{file}.fa"
    output:
        temp("{file}_clean.fa")
    shell:
        "sed 's/_[^_]*$//' {input} > {output}"


rule isoformswitchanalyser_importData_flair:
    input:
        counts = "flair/{dataset}/flair_{dataset}_counts_matrix.tsv",
        sample_info = "sample_info/sampleInfo.csv",
        gtf = "flair/{dataset}/flair.collapse.isoforms.gtf",
        transcripts = "flair/{dataset}/flair.collapse.isoforms_clean.fa",
        test = "res/drimseq/{dataset}_flair/{dataset}_flair_drimSeqStageR.csv"
    output:
        "res/drimseq/{dataset}_flair/{dataset}_flair_sal.rds"
    wildcard_constraints:
        dataset = "cdna|teloprime"
    script:
        "drimseq_isoformswitchanalyser_importData_flair.R"


def get_txdb(wildcards):
    if wildcards.dataset in ["teloprime", "illumina", "cdna"]:
        txdb = "annotation/annotation_txdb.sqlite"
    elif wildcards.dataset == "teloprime_flair":
        txdb = "flair/teloprime/flair.collapse.isoforms_txdb.sqlite"
    elif wildcards.dataset == "cdna_flair":
        txdb = "flair/cdna/flair.collapse.isoforms_txdb.sqlite"
    return txdb


def get_axis(wildcards):
    if wildcards.dataset == "illumina":
        axis = "illumina"
    else:
        axis = "ont"
    return axis


rule drimseq_report:
    input:
        dmds = "res/drimseq/{dataset}/{file}_dmds_dtu.rds",
        res = "res/drimseq/{dataset}/{file}_drimSeqStageR.csv",
        biomaRt_gene = "annotation/biomaRt_gene.rds",
        biomaRt_tx = "annotation/biomaRt_tx.rds",
        switchList = "res/drimseq/{dataset}/{dataset}_sal.rds"
    output:
        "res/drimseq/{dataset}/{file}_drimSeqStageR.html"
    params:
        axis = get_axis
    wildcard_constraints:
        dataset = "teloprime|cdna|illumina|cdna_flair|teloprime_flair"
    script:
        "drimseq_stagerAnalysis.Rmd"
