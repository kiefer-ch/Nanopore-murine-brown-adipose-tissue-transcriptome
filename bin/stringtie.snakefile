rule stringtie:
    input:
        bam = "bam/illumina/{sample}_Aligned.sortedByCoord.out.bam",
        annotation = "annotation/annotation.gtf"
    output:
        "stringtie/illumina/stringtie_illumina_{sample}.gtf"
    threads:
        4
    params:
        label = "{sample}"
    shell:
        "stringtie \
            {input.bam} \
            --rf \
            -p {threads} \
            -l {params.label} \
            -G {input.annotation} \
            -o {output} \
            -j 10"


rule stringtie_merge:
    input:
        gtfs = expand("stringtie/illumina/stringtie_illumina_{sample}.gtf",
            sample=SAMPLES_ont)
    output:
        "stringtie/illumina/stringtie_illumina_merged.gtf"
    params:
    threads:
        8
    shell:
        "stringtie --merge \
            -p {threads} \
            -l 'stringtie_merge' \
            -o {output} \
            {input.gtfs}"


rule stringtie_gffcompare:
    input:
        "annotation/genome.fa.fai",
        stringtie = "stringtie/illumina/stringtie_illumina_merged.gtf",
        reference = "annotation/annotation.gtf",
        genome = "annotation/genome.fa"
    output:
        "stringtie/illumina/gffcompare/illumina_gffcompare.loci",
        "stringtie/illumina/gffcompare/illumina_gffcompare.stats",
        "stringtie/illumina/gffcompare/illumina_gffcompare.annotated.gtf",
        "stringtie/illumina/gffcompare/illumina_gffcompare.tracking"
    params:
        out_prefix = "stringtie/illumina/gffcompare/illumina_gffcompare"
    shell:
        "gffcompare -R -T \
            -o  {params.out_prefix} \
            -r {input.reference} \
            {input.stringtie}"


rule stringtie_sqanti:
    input:
        "annotation/genome.fa.fai",
        isoforms = "stringtie/illumina/stringtie_illumina_merged.gtf",
        annotation = "annotation/annotation.gtf",
        genome = "annotation/genome.fa"
    output:
        "stringtie/illumina/sqanti/stringtie.collapse.isoforms_report.pdf",
        "stringtie/illumina/sqanti/stringtie.collapse.isoforms_classification.txt"
    params:
        out_dir = "stringtie/illumina/sqanti"
    threads:
        10
    shell:
        "sqanti_qc2 -g \
            -d {params.out_dir} \
            {input.isoforms} {input.annotation} {input.genome}"
