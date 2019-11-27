# deseq
rule tximport_deseq_gene:
    input:
        salmon_out = expand("salmon/{sample}/quant.sf", sample=SAMPLES),
        txdb = "annotation/annotation_txdb.sqlite",
        sample_info = "sample_info/sampleInfo.csv"
    params:
        txOut = 0,
        design = "~condition_temp"
    output:
        "res/deseq/illumina/genelevel_all/illumina_genelevel_all_dds.rds"
    script:
        "txImport_deseq.R"

rule tximport_deseq_transcript:
    input:
        salmon_out = expand("salmon/{sample}/quant.sf", sample=SAMPLES),
        txdb = "annotation/annotation_txdb.sqlite",
        sample_info = "sample_info/sampleInfo.csv"
    params:
        txOut = 1,
        design = "~condition_temp"
    output:
        "res/deseq/illumina/txlevel_all/illumina_txlevel_all_dds.rds"
    script:
        "txImport_deseq.R"

# dexseq
rule tximport_dexseq_illumina:
    input:
        salmone_out = expand("salmon/{sample}/quant.sf", sample=SAMPLES_ont),
        sample_info = "sample_info/sampleInfo.csv"
    output:
        "res/dexseq/illumina/dexseq_scaledTPM.rds"
    script:
        "txImport_dexseq.R"
