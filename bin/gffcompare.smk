# gffcompare
rule gffcompare:
    input:
        gtfs =
            [expand("data/reannotation/flair/annotation/{dataset}_flair.isoforms.gtf",
                    dataset=["cdna", "teloprime", "rna"]),
                expand("data/reannotation/stringtie/{dataset}_stringtie.gtf",
                    dataset=["cdna", "teloprime", "rna", "illumina"])],
        reference = "data/annotation/annotation.gtf"
    output:
        multiext("data/comparisons/reannotation/gffcompare/R/gffcmp", ".loci",
            ".stats", ".tracking"),
        "data/comparisons/reannotation/gffcompare/R/gffcmp.combined.gtf",
        expand("data/reannotation/stringtie/gffcmp.{dataset}_stringtie.gtf.tmap",
            dataset=["cdna", "teloprime", "rna", "illumina"]),
        expand("data/reannotation/flair/annotation/gffcmp.{dataset}_flair.isoforms.gtf.tmap",
               dataset=["cdna", "teloprime", "rna"])
    shell:
        """
        gffcompare \
            -r {input.reference} -R \
            {input.gtfs} &&
        mv gffcmp.* data/comparisons/reannotation/gffcompare/R
        """


rule gffcompare_getPrecision:
# -Q flag is set to ignore novel transcripts
# -T is set so no tmap and refmap files are produced
    input:
        gtfs =
            [expand("data/reannotation/flair/annotation/{dataset}_flair.isoforms.gtf",
                    dataset=["cdna", "teloprime", "rna"]),
                expand("data/reannotation/stringtie/{dataset}_stringtie.gtf",
                    dataset=["cdna", "teloprime", "rna", "illumina"])],
        reference = "data/annotation/annotation.gtf"
    output:
        multiext("data/comparisons/reannotation/gffcompare/RQ/gffcmp", ".loci",
            ".stats", ".tracking"),
        "data/comparisons/reannotation/gffcompare/RQ/gffcmp.combined.gtf"
    shell:
        """
        gffcompare \
            -r {input.reference} -R -Q -T \
            {input.gtfs} &&
        mv gffcmp.* data/comparisons/reannotation/gffcompare/RQ
        """


# squanti
rule filterUnknownStrand_fasta:
    input:
        "{file}.gtf"
    output:
        "{file}_noUnknownStrand.gtf"
    conda:
        "../envs/r_4.1.2.yaml"
    script:
        "gffcompare_filterFasta.R"


rule sqanti_stringtie:
    input:
        "data/annotation/genome.fa.fai",
        isoforms = "data/reannotation/stringtie/{dataset}_stringtie_noUnknownStrand.gtf",
        annotation = "data/annotation/annotation.gtf",
        genome = "data/annotation/genome.fa"
    output:
        "data/comparisons/reannotation/squanti/stringtie/{dataset}/{dataset}_stringtie_noUnknownStrand_sqanti_report.pdf",
        "data/comparisons/reannotation/squanti/stringtie/{dataset}/{dataset}_stringtie_noUnknownStrand_classification.txt"
    params:
        out_dir = "data/comparisons/reannotation/squanti/stringtie/{dataset}"
    shell:
        "sqanti_qc2 -g \
            -d {params.out_dir} \
            {input.isoforms} {input.annotation} {input.genome}"


rule sqanti_flair:
    input:
        "data/annotation/genome.fa.fai",
        isoforms = "data/reannotation/flair/annotation/{dataset}_flair.isoforms.gtf",
        annotation = "data/annotation/annotation.gtf",
        genome = "data/annotation/genome.fa"
    output:
        "data/comparisons/reannotation/squanti/flair/{dataset}/{dataset}_flair.isoforms_sqanti_report.pdf",
        "data/comparisons/reannotation/squanti/flair/{dataset}/{dataset}_flair.isoforms_classification.txt"
    params:
        out_dir = "data/comparisons/reannotation/squanti/flair/{dataset}"
    shell:
        "sqanti_qc2 -g \
            -d {params.out_dir} \
            {input.isoforms} {input.annotation} {input.genome}"
