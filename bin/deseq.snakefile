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


rule import_deseq_cdna_gene:
    input:
        counts = expand("res/wien/direct_cDNA/pool{pool}/20200108_pool{pool}_transcriptome_barcode{barcode}_q7_counts.tsv.gz",
            barcode=["07","08","09","10","11","12"], pool=[1,2]),
        sample_info = "sample_info/sampleInfo.csv"
    params:
        txOut = 0,
        design = "~condition_temp"
    output:
        dds = "res/deseq/cdna/genelevel/cdna_genelevel_dds.rds"
    script:
        "deseq_cDNA.R"


rule import_deseq_cdna_tx:
    input:
        counts = expand("res/wien/direct_cDNA/pool{pool}/20200108_pool{pool}_transcriptome_barcode{barcode}_q7_counts.tsv.gz",
            barcode=["07","08","09","10","11","12"], pool=[1,2]),
        sample_info = "sample_info/sampleInfo.csv"
    params:
        txOut = 1,
        design = "~condition_temp"
    output:
        dds = "res/deseq/cdna/txlevel/cdna_txlevel_dds.rds"
    script:
        "deseq_cDNA.R"


rule import_deseq_rna_gene:
    input:
        counts = expand("res/wien/direct_RNA/ChrKiefer_20200306/counts_per_transcript/transcriptome_{temp}_q7_counts.tsv.gz",
            temp=["cool", "rt"])
    params:
        txOut = 0,
        design = "~condition_temp"
    output:
        dds = "res/deseq/rna/genelevel/rna_genelevel_dds.rds"
    script:
        "deseq_rna.R"


rule import_deseq_rna_tx:
    input:
        counts = expand("res/wien/direct_RNA/ChrKiefer_20200306/counts_per_transcript/transcriptome_{temp}_q7_counts.tsv.gz",
            temp=["cool", "rt"])
    params:
        txOut = 1,
        design = "~condition_temp"
    output:
        dds = "res/deseq/rna/txlevel/rna_txlevel_dds.rds"
    script:
        "deseq_rna.R"


# qc
rule deseq_qc:
    input:
        dds = "{method}_dds.rds"
    output:
        "{method}_qc.html"
    script:
        "deseq_qc.Rmd"


# cm export
def get_biomart(wildcards):
    if wildcards.level == "genelevel":
        return "annotation/biomaRt_gene.rds"
    elif wildcards.level == "txlevel":
        return "annotation/biomaRt_tx.rds"


rule deseq_exportCm_illumina:
    input:
        dds = "{file_path}/illumina_{level}_dds.rds",
        biomart = get_biomart
    params:
        tpm = 1,
        level = "{level}",
        vst = 0
    output:
        cts = "{file_path}/illumina_{level}_cm_cts.csv.gz",
        ntd = "{file_path}/illumina_{level}_cm_ntd.csv.gz",
        rld = "{file_path}/illumina_{level}_cm_rld.csv.gz",
        tpm = "{file_path}/illumina_{level}_cm_tpm.csv.gz"
    wildcard_constraints:
        level = "genelevel|txlevel",
    script:
        "deseq_exportCm.R"


rule deseq_exportCm_ont:
    input:
        dds = "{file_path}/{dataset}_{level}_dds.rds",
        biomart = get_biomart
    params:
        tpm = 0,
        level = "{level}",
        vst = 0
    wildcard_constraints:
        dataset = "teloprime|cdna",
        level = "genelevel|txlevel"
    output:
        cts = "{file_path}/{dataset}_{level}_cm_cts.csv.gz",
        ntd = "{file_path}/{dataset}_{level}_cm_ntd.csv.gz",
        rld = "{file_path}/{dataset}_{level}_cm_rld.csv.gz"
    script:
        "deseq_exportCm.R"

# DGE
def get_threshold(wildcards):
    if wildcards.shrink == "apeglm":
        return 1
    elif wildcards.shrink == "noShrink":
        return 0

rule deseq_dge:
    input:
        dds = "{file_path}/{dataset}_{level}_dds.rds",
        biomart = get_biomart
    output:
        "{file_path}/{dataset}_{level}_{shrink}_results.csv.gz",
    params:
        level = "{level}",
        shrink = "{shrink}",
        lfcThreshold = get_threshold,
        alpha = 0.05
    wildcard_constraints:
        level = "genelevel|txlevel",
        shrink = "apeglm|noShrink"
    threads: 4
    script:
        "deseq_dge.R"


# pathway analysis
rule topgo_analysis:
    input:
        "{file}_{level}_{shrink}_results.csv.gz"
    output:
        "{file}_{level}_{shrink}_topgo.rds"
    params:
        level = "{level}",
        shrink = "{shrink}"
    wildcard_constraints:
        level = "genelevel|txlevel",
        shrink = "apeglm|noShrink"
    script:
        "deseq_topgo.R"


rule reactome_analysis:
    input:
        "{file}_{level}_{shrink}_results.csv.gz"
    output:
        "{file}_{level}_{shrink}_reactome.rds",
        "{file}_{level}_{shrink}_reactomeGL.rds"
    params:
        level = "{level}",
        shrink = "{shrink}"
    wildcard_constraints:
        level = "genelevel|txlevel",
        shrink = "apeglm|noShrink"
    script:
        "deseq_reactome.R"


# report
rule deseq_report:
    input:
        dds = "{file}_{level}_dds.rds",
        rld = "{file}_{level}_cm_rld.csv.gz",
        res = "{file}_{level}_{shrink}_results.csv.gz",
        topgo = "{file}_{level}_{shrink}_topgo.rds",
        reactome = "{file}_{level}_{shrink}_reactome.rds",
        reactome_genelist = "{file}_{level}_{shrink}_reactomeGL.rds",
        grouped_biotypes = "data/biotype_groups.csv"
    params:
        level = "{level}",
        shrink = "{shrink}",
        lfcThreshold = get_threshold,
        alpha = .05
    wildcard_constraints:
        level = "genelevel|txlevel",
        shrink = "apeglm|noShrink"
    output:
        "{file}_{level}_{shrink}_report.html"
    script:
        "deseq_report.Rmd"
