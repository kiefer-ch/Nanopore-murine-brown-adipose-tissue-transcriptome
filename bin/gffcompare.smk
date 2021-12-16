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
        multiext("data/comparisons/reannotation/gffcompare/gffcompare", ".loci",
            ".stats", ".tracking"),
        "data/comparisons/reannotation/gffcompare/gffcompare.combined.gtf"
#        expand("data/reannotation/stringtie/gffcompare.{dataset}_stringtie.gtf.tmap",
#            dataset=["illumina", "cdna", "teloprime", "rna"]),
#        expand("data/reannotation/flair/annotation/gffcompare.{dataset}_flair.isoforms.gtf.tmap",
#               dataset=["cdna", "teloprime", "rna"])
    params:
        out_prefix = "data/comparisons/reannotation/gffcompare/gffcompare"
    shell:
        "gffcompare -R \
            -o {params.out_prefix} \
            -r {input.reference} \
            {input.gtfs}"


# squanti
rule filterUnknownStrand_fasta:
    input:
        "{file}.gtf"
    output:
        temp("{file}_noUnknownStrand.gtf")
    script:
        "gffcompare_filterFasta.R"


rule sqanti_stringtie:
    input:
        "data/annotation/genome.fa.fai",
        isoforms = "data/reannotation/stringtie/{dataset}_stringtie_noUnknownStrand.gtf",
        annotation = "data/annotation/annotation.gtf",
        genome = "data/annotation/genome.fa"
    output:
        "data/comparisons/reannotation/squanti/{dataset}/{dataset}_stringtie_noUnknownStrand_sqanti_report.pdf",
        "data/comparisons/reannotation/squanti/{dataset}/{dataset}_stringtie_noUnknownStrand_classification.txt"
    params:
        out_dir = "data/comparisons/reannotation/squanti/{dataset}"
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
        "data/comparisons/reannotation/squanti/{dataset}/{dataset}_flair.isoforms_sqanti_report.pdf",
        "data/comparisons/reannotation/squanti/{dataset}/{dataset}_flair.isoforms_classification.txt"
    params:
        out_dir = "data/comparisons/reannotation/squanti/{dataset}"
    shell:
        "sqanti_qc2 -g \
            -d {params.out_dir} \
            {input.isoforms} {input.annotation} {input.genome}"
