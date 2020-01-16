# data import
rule tximport_deseq_illumina_gene:
    input:
        salmon_out = expand("salmon/{sample}/quant.sf", sample=SAMPLES),
        txdb = "annotation/annotation_txdb.sqlite",
        sample_info = "sample_info/sampleInfo.csv"
    params:
        txOut = 0,
        design = "~condition_temp"
    output:
        "res/deseq/illumina/genelevel/illumina_genelevel_dds.rds"
    script:
        "deseq_txImport.R"


rule tximport_deseq_illumina_transcript:
    input:
        salmon_out = expand("salmon/{sample}/quant.sf", sample=SAMPLES),
        txdb = "annotation/annotation_txdb.sqlite",
        sample_info = "sample_info/sampleInfo.csv"
    params:
        txOut = 1,
        design = "~condition_temp"
    output:
        "res/deseq/illumina/txlevel/illumina_txlevel_dds.rds"
    script:
        "deseq_txImport.R"


rule import_deseq_teloprime_gene:
    input:
        counts = "res/wien/teloprime/DE/ONT_newbasecalling/6samples_counts.tsv",
        sample_info = "sample_info/sampleInfo.csv"
    params:
        txOut = 0,
        design = "~condition_temp"
    output:
        dds = "res/deseq/teloprime/genelevel/teloprime_genelevel_dds.rds"
    script:
        "deseq_teloprime.R"


rule import_deseq_teloprime_transcript:
    input:
        counts = "res/wien/teloprime/DE/ONT_newbasecalling/6samples_counts.tsv",
        sample_info = "sample_info/sampleInfo.csv"
    params:
        txOut = 1,
        design = "~condition_temp"
    output:
        dds = "res/deseq/teloprime/txlevel/teloprime_txlevel_dds.rds"
    script:
        "deseq_teloprime.R"


# qc
rule deseq_qc:
    input:
        dds = "{method}_dds.rds"
    output:
        "{method}_qc.html"
    script:
        "deseq_qc.Rmd"


# cm export
rule deseq_exportCm_illumina_genelevel:
    input:
        dds = "{file_path}/illumina_genelevel_dds.rds",
        biomart = "annotation/biomaRt_gene.rds"
    params:
        tpm = 1,
        level = "genelevel",
        vst = 0
    output:
        cts = "{file_path}/illumina_genelevel_cm_cts.csv.gz",
        ntd = "{file_path}/illumina_genelevel_cm_ntd.csv.gz",
        rld = "{file_path}/illumina_genelevel_cm_rld.csv.gz",
        tpm = "{file_path}/illumina_genelevel_cm_tpm.csv.gz"
    script:
        "deseq_exportCm.R"


rule deseq_exportCm_illumina_txlevel:
    input:
        dds = "{file_path}/illumina_txlevel_dds.rds",
        biomart = "annotation/biomaRt_tx.rds"
    params:
        tpm = 1,
        level = "txlevel",
        vst = 0
    output:
        cts = "{file_path}/illumina_txlevel_cm_cts.csv.gz",
        ntd = "{file_path}/illumina_txlevel_cm_ntd.csv.gz",
        rld = "{file_path}/illumina_txlevel_cm_rld.csv.gz",
        tpm = "{file_path}/illumina_txlevel_cm_tpm.csv.gz"
    script:
        "deseq_exportCm.R"

rule deseq_exportCm_ont_genelevel:
    input:
        dds = "{file_path}/{dataset}_genelevel_dds.rds",
        biomart = "annotation/biomaRt_gene.rds"
    params:
        tpm = 0,
        level = "genelevel",
        vst = 0
    wildcard_constraints:
        datatset = "teloprime"
    output:
        cts = "{file_path}/{dataset}_genelevel_cm_cts.csv.gz",
        ntd = "{file_path}/{dataset}_genelevel_cm_ntd.csv.gz",
        rld = "{file_path}/{dataset}_genelevel_cm_rld.csv.gz"
    script:
        "deseq_exportCm.R"


rule deseq_exportCm_ont_txlevel:
    input:
        dds = "{file_path}/{dataset}_txlevel_dds.rds",
        biomart = "annotation/biomaRt_tx.rds"
    params:
        tpm = 0,
        level = "txlevel",
        vst = 0
    output:
        cts = "{file_path}/{dataset}_txlevel_cm_cts.csv.gz",
        ntd = "{file_path}/{dataset}_txlevel_cm_ntd.csv.gz",
        rld = "{file_path}/{dataset}_txlevel_cm_rld.csv.gz"
    script:
        "deseq_exportCm.R"


# differential expression analysis


# DGE
rule deseq_genelevel_ont:
    threads: 4
    params:
        out_folder = "res/deseq/illumina/genelevel_ont"
    input:
        "data/dds_gencode.vM22_gene_ont.rds",
        "data/biotype_groups.csv"
    output:
        "res/deseq/illumina/genelevel_ont/deseq_genelevel_ont.html"
    script:
        "deseq_genelevel_ont.Rmd"

rule deseq_genelevel_all:
    threads: 4
    params:
        out_folder = "res/deseq/illumina/genelevel_all"
    input:
        "data/dds_gencode.vM22_gene.rds",
        "data/biotype_groups.csv"
    output:
        "res/deseq/illumina/genelevel_all/deseq_genelevel_all.html"
    script:
        "deseq_genelevel_all.Rmd"

# DTE
rule deseq_txlevel_ont:
    threads: 4
    params:
        out_folder = "res/deseq/illumina/txlevel_ont"
    input:
        "data/dds_gencode.vM22_transcript_ont.rds",
        "data/biotype_groups.csv"
    output:
        "res/deseq/illumina/txlevel_ont/deseq_txlevel_ont.html"
    script:
        "deseq_txlevel_ont.Rmd"

rule deseq_txlevel_all:
    threads: 4
    params:
        out_folder = "res/deseq/illumina/txlevel_all"
    input:
        "data/dds_gencode.vM22_transcript.rds",
        "data/biotype_groups.csv"
    output:
        "res/deseq/illumina/txlevel_all/deseq_txlevel_all.html"
    script:
        "deseq_txlevel_all.Rmd"
