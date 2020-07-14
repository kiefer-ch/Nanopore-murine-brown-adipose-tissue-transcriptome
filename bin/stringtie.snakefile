STRINGTIE = config["STRINGTIE"]

rule stringtie_illumina:
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
        """
        {STRINGTIE} \
            {input.bam} \
            --rf \
            -p {threads} \
            -l {params.label} \
            -G {input.annotation} \
            -o {output} \
            -j 10
        """


rule stringtie_ont:
    input:
        bam = "bam/{dataset}/{sample}_genome.bam",
        annotation = "annotation/annotation.gtf"
    output:
        "stringtie/{dataset}/stringtie_{dataset}_{sample}.gtf"
    threads:
        4
    params:
        label = "{sample}"
    shell:
        """
        {STRINGTIE} \
            {input.bam} \
            -L \
            -p {threads} \
            -l {params.label} \
            -G {input.annotation} \
            -o {output} \
            -j 10
        """


def get_stringtie_merge_gtf_names(wildcards):
    if wildcards.dataset == "illumina":
        gtfs = expand("stringtie/illumina/stringtie_illumina_{sample}.gtf",
            sample=SAMPLE_INFO_ont["illumina"].tolist())
    elif wildcards.dataset == "cdna":
        gtfs = expand("stringtie/cdna/stringtie_cdna_{sample}.gtf",
            sample=SAMPLE_INFO_ont["cdna"].tolist())
    elif wildcards.dataset == "teloprime":
        gtfs = expand("stringtie/teloprime/stringtie_teloprime_{sample}.gtf",
            sample=SAMPLE_INFO_ont["ont"].tolist())
    return gtfs


rule stringtie_merge:
    input:
        gtfs = get_stringtie_merge_gtf_names
    output:
        "stringtie/{dataset}/stringtie_{dataset}_merged.gtf"
    threads:
        8
    params:
        label = "stringtie_merge_{dataset}"
    shell:
        """
        {STRINGTIE} --merge \
            -p {threads} \
            -l {params.label} \
            -o {output} \
            {input.gtfs}
        """
