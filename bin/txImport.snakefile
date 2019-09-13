rule tximport_gene:
    input:
        expand("salmon/{sample}/quant.sf", sample=SAMPLES),
        annotation = "annotation/annotation.gtf",
        sample_info = "sample_info/sampleInfo.csv"
    output:
        "data/dds_gencode.vM22_gene.rds"
    shell:
        "bin/txImport.R {input.sample_info} {input.annotation} {output} --genelevel"

rule tximport_transcript:
    input:
        expand("salmon/{sample}/quant.sf", sample=SAMPLES),
        annotation = "annotation/annotation.gtf",
        sample_info = "sample_info/sampleInfo.csv"
    output:
        "data/dds_gencode.vM22_transcript.rds"
    shell:
        "bin/txImport.R {input.sample_info} {input.annotation} {output} --txlevel"

rule tximport_all:
    input:
        "data/dds_gencode.vM22_transcript.rds",
        "data/dds_gencode.vM22_gene.rds"

