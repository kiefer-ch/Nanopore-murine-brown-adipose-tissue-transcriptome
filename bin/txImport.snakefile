rule tximport_gene:
    input:
        expand("salmon/{sample}/quant.sf", sample=SAMPLES),
        annotation = "annotation/annotation.gtf",
        sample_info = "sample_info/sampleInfo.csv"
    output:
        "data/dds_gencode.vM22_gene.rds"
    shell:
        "bin/txImport.R {input.sample_info} {input.annotation} {output} --genelevel illumina"

rule tximport_transcript:
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

rule tximport_dtu:
    input:
        expand("salmon/{sample}/quant.sf", sample=SAMPLES_ont),
        "sample_info/sampleInfo.csv"
    output:
        "data/scaledTPM.rds"
    shell:
        "bin/txImport_DTU.R"

rule tximport_all:
    input:
        "data/dds_gencode.vM22_transcript.rds",
        "data/dds_gencode.vM22_gene.rds",
        "data/dds_gencode.vM22_transcript_ont.rds",
        "data/dds_gencode.vM22_gene_ont.rds"
