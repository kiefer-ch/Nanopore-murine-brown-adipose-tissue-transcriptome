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
    conda:
        "../envs/stringtie_2.2.0.yaml"
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
    conda:
        "../envs/stringtie_2.2.0.yaml"
    shell:
        """
        gffcompare \
            -r {input.reference} -R -Q -T \
            {input.gtfs} &&
        mv gffcmp.* data/comparisons/reannotation/gffcompare/RQ
        """