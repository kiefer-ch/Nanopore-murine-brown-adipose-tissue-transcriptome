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
