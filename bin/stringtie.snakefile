rule stringtie:
    input:
        bam = "bam/illumina/{sample}_Aligned.sortedByCoord.out.bam",
        annotation = "annotation/annotation.gtf"
    output:
        "stringtie/illumina//illumina/stringtie_illumina_{sample}.gtf"
    threads:
        4
    shell:
        "stringtie \
            {input.bam} \
            --rf \
            -p 4 \
            -l {sample} \
            -G {input.annotation} \
            -o {output} \
            -j 5"
