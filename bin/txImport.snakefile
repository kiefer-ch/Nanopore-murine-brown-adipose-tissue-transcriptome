rule tximport_salmon_gene:
    input:
        salmon_out = expand("salmon/{sample}/quant.sf", sample=SAMPLES),
        txdb = "annotation/annotation_txdb.sqlite",
        sample_info = "sample_info/sampleInfo.csv"
    params:
        level = "gene"
    output:
        "res/deseq/illumina/genelevel_all/illumina_genelevel_all_dds.rds"
    shell:
        "bin/txImport.R"

rule tximport_salmon_transcript:
    input:
        expand("salmon/{sample}/quant.sf", sample=SAMPLES),
        annotation = "annotation/annotation.gtf",
        sample_info = "sample_info/sampleInfo.csv"
    output:
        "data/dds_gencode.vM22_transcript.rds"
    shell:
        "bin/txImport.R {input.sample_info} {input.annotation} {output} --txlevel illumina"

rule tximport_dtu_ont:
    input:
        expand("salmon/{sample}/quant.sf", sample=SAMPLES_ont),
        "sample_info/sampleInfo.csv"
    output:
        "data/scaledTPM_ont.rds"
    shell:
        "bin/txImport_DTU.R {output} ont"

rule tximport_dtu_all:
    input:
        expand("salmon/{sample}/quant.sf", sample=SAMPLES_ont),
        "sample_info/sampleInfo.csv"
    output:
        "data/scaledTPM_all.rds"
    shell:
        "bin/txImport_DTU.R {output} illumina"
