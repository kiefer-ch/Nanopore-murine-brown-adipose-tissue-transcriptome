rule stringtie_illumina:
    input:
        bam = "data/bam/illumina/{sample}_Aligned.sortedByCoord.out.bam",
        annotation = "data/annotation/annotation.gtf"
    output:
        "data/reannotation/stringtie/illumina/illumina_{sample}_stringtie.gtf"
    threads:
        4
    params:
        label = "{sample}",
        SJ_cutoff = config["SJ_cutoff"]
    conda:
        "../envs/stringtie_2.2.0.yaml"
    shell:
        """
        stringtie \
            {input.bam} \
            --rf \
            -p {threads} \
            -l {params.label} \
            -G {input.annotation} \
            -o {output} \
            -j {params.SJ_cutoff} \
            -c 1 \
            -f 0.01
        """


def get_illumina_bam(wildcards):
    if wildcards.dataset == "cdna":
        illumina = SAMPLE_INFO_ont.set_index("cdna").loc[wildcards.barcode, "illumina"]
        filename = "data/bam/illumina/{}_Aligned.sortedByCoord.out.bam".format(illumina)
    elif wildcards.dataset == "teloprime":
        illumina = SAMPLE_INFO_ont.set_index("ont").loc[wildcards.barcode, "illumina"]
        filename = "data/bam/illumina/{}_Aligned.sortedByCoord.out.bam".format(illumina)
    elif wildcards.dataset == "rna":
        if wildcards.barcode == "rt":
            filename = "data/bam/illumina/5034_S33_Aligned.sortedByCoord.out.bam"
        if wildcards.barcode == "cool":
            filename = "data/bam/illumina/5035_S34_Aligned.sortedByCoord.out.bam"
    return filename


rule stringtie_ont:
    input:
        bam_long = "data/bam/{dataset}/merged/{dataset}_{barcode}_genome_primaryOnly.bam",
        bam_short = get_illumina_bam,
        annotation = "data/annotation/annotation.gtf"
    output:
        "data/reannotation/stringtie/{dataset}/{dataset}_{barcode}_stringtie.gtf"
    threads:
        4
    params:
        label = "{barcode}",
        SJ_cutoff = config["SJ_cutoff"]
    wildcard_constraints:
        dataset = "teloprime|cdna|rna"
    conda:
        "../envs/stringtie_2.2.0.yaml"
    shell:
        """
        stringtie \
            {input.bam_short} {input.bam_long} \
            --mix \
            -p {threads} \
            -l {params.label} \
            -G {input.annotation} \
            -o {output} \
            -j {params.SJ_cutoff} \
            -c 1 \
            -f 0.01
        """


def get_stringtie_merge_gtf_names(wildcards):
    if wildcards.dataset == "illumina":
        gtfs = expand("data/reannotation/stringtie/illumina/illumina_{sample}_stringtie.gtf",
            sample=SAMPLE_INFO_ont["illumina"].tolist())
    elif wildcards.dataset == "cdna":
        gtfs = expand("data/reannotation/stringtie/cdna/cdna_{sample}_stringtie.gtf",
            sample=SAMPLE_INFO_ont["cdna"].tolist())
    elif wildcards.dataset == "teloprime":
        gtfs = expand("data/reannotation/stringtie/teloprime/teloprime_{sample}_stringtie.gtf",
            sample=SAMPLE_INFO_ont["ont"].tolist())
    elif wildcards.dataset == "rna":
        gtfs = expand("data/reannotation/stringtie/rna/rna_{sample}_stringtie.gtf",
            sample=["cool", "rt"])
    return gtfs


rule stringtie_merge:
    input:
        gtfs = get_stringtie_merge_gtf_names
    output:
        "data/reannotation/stringtie/{dataset}_stringtie.gtf"
    threads:
        8
    params:
        label = "stringtie_merge_{dataset}",
        SUPPORT_cutoff = config["SUPPORT_cutoff"]
    wildcard_constraints:
        dataset = "teloprime|cdna|rna|illumina"
    conda:
        "../envs/stringtie_2.2.0.yaml"
    shell:
        """
        stringtie --merge \
            -p {threads} \
            -l {params.label} \
            -c {params.SUPPORT_cutoff} \
            -f 0.05 \
            -o {output} \
            {input.gtfs}
        """


rule stringtie_getFasta:
    input:
        "data/annotation/genome.fa.fai",
        gtf = "data/reannotation/stringtie/{dataset}_stringtie_noUnknownStrand.gtf",
        genome = "data/annotation/genome.fa"
    output:
        "data/reannotation/stringtie/{dataset}_stringtie_noUnknownStrand.fa"
    wildcard_constraints:
        dataset = "illumina|cdna|teloprime"
    conda:
        "../envs/stringtie_2.2.0.yaml"
    shell:
        "gffread --ignore-locus -w {output} -g {input.genome} {input.gtf}"
