rule tximport_gene_illumina:
    input:
        expand("salmon/{sample}/quant.sf", sample=SAMPLES),
        annotation = "annotation/annotation.gtf",
        sample_info = "sample_info/sampleInfo.csv"
    output:
        "data/dds_gencode.vM22_gene.rds"
    shell:
        "bin/txImport.R {input.sample_info} {input.annotation} {output} --genelevel illumina"

rule tximport_transcript_illumina:
    input:
        expand("salmon/{sample}/quant.sf", sample=SAMPLES),
        annotation = "annotation/annotation.gtf",
        sample_info = "sample_info/sampleInfo.csv"
    output:
        "data/dds_gencode.vM22_transcript.rds"
    shell:
        "bin/txImport.R {input.sample_info} {input.annotation} {output} --txlevel illumina"

rule tximport_gene_ont:
    input:
        expand("salmon/{sample}/quant.sf", sample=SAMPLES_ont),
        annotation = "annotation/annotation.gtf",
        sample_info = "sample_info/sampleInfo.csv"
    output:
        "data/dds_gencode.vM22_gene_ont.rds"
    shell:
        "bin/txImport.R {input.sample_info} {input.annotation} {output} --genelevel ont"

rule tximport_transcript_ont:
    input:
        expand("salmon/{sample}/quant.sf", sample=SAMPLES_ont),
        annotation = "annotation/annotation.gtf",
        sample_info = "sample_info/sampleInfo.csv"
    output:
        "data/dds_gencode.vM22_transcript_ont.rds"
    shell:
        "bin/txImport.R {input.sample_info} {input.annotation} {output} --txlevel ont"

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
