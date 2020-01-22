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
rule deseq_dge_apeglm:
    input:
        dds = "{file_path}/{dataset}_genelevel_dds.rds",
        biomart = "annotation/biomaRt_gene.rds"
    output:
        "{file_path}/{dataset}_genelevel_apeglm_results.csv.gz",
    params:
        level = "genelevel",
        shrink = 1,
        lfcThreshold = 1,
        alpha = 0.05
    threads: 4
    script:
        "deseq_dge.R"


rule deseq_dge_noShrink:
    input:
        dds = "{file_path}/{dataset}_genelevel_dds.rds",
        biomart = "annotation/biomaRt_gene.rds"
    output:
        "{file_path}/{dataset}_genelevel_noshrink_results.csv.gz",
    params:
        level = "genelevel",
        shrink = 0,
        lfcThreshold = 1,
        alpha = 0.05
    threads: 4
    script:
        "deseq_dge.R"


# DTE
rule deseq_dte_apeglm:
    input:
        dds = "{file_path}/{dataset}_txlevel_dds.rds",
        biomart = "annotation/biomaRt_tx.rds"
    output:
        "{file_path}/{dataset}_txlevel_apeglm_results.csv.gz",
    params:
        level = "txlevel",
        shrink = 1,
        lfcThreshold = 1,
        alpha = 0.05
    threads: 4
    script:
        "deseq_dge.R"


rule deseq_dte_noShrink:
    input:
        dds = "{file_path}/{dataset}_txlevel_dds.rds",
        biomart = "annotation/biomaRt_tx.rds"
    output:
        "{file_path}/{dataset}_txlevel_noshrink_results.csv.gz",
    params:
        level = "txlevel",
        shrink = 0,
        lfcThreshold = 1,
        alpha = 0.05
    threads: 4
    script:
        "deseq_dge.R"


# pathway analysis
rule topgo_analysis:
    input:
        "{file}_results.csv.gz"
    output:
        "{file}_topgo.csv.gz"
    params:
        level = "genelevel",
        shrink = 1
    script:
        "deseq_topgo.R"


rule reactome_analysis:
    input:
        "{file}_{level}_{shrink}_results.csv.gz"
    output:
        "{file}_{level}_{shrink}_reactome.rds"
    params:
        level = "{level}",
        shrink = "{shrink}"
    wildcard_constraints:
        level = "genelevel|txlevel",
        shrink = "apeglm|noShrink"
    script:
        "deseq_reactome.R"


# summary
rule deseq_summary:
    input:
        res = "{file}_results.csv.gz",
        topgo = "{file}_topgo.csv.gz",
        reactome = "{file}_reactome.rds"
    output:
        "{file}.html"
    script:
        "deseq_summary.Rmd"
