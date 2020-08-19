rule gffcompare:
    input:
        gtfs = ["stringtie/illumina/stringtie_illumina_merged.gtf",
                expand("flair/{dataset}/flair.collapse.{dataset}.isoforms.gtf",
                       dataset=["cdna", "teloprime", "rna"])],
        reference = "annotation/annotation.gtf"
    output:
        "res/comparisons/gffcompare/all_gffcompare.loci",
        "res/comparisons/gffcompare/all_gffcompare.stats",
        "res/comparisons/gffcompare/all_gffcompare.combined.gtf",
        "res/comparisons/gffcompare/all_gffcompare.tracking",
        "stringtie/illumina/all_gffcompare.stringtie_illumina_merged.gtf.tmap",
        expand("flair/{dataset}/all_gffcompare.flair.collapse.{dataset}.isoforms.gtf.tmap",
               dataset=["cdna", "teloprime", "rna"])
    params:
        out_prefix = "res/comparisons/gffcompare/all_gffcompare"
    shell:
        "gffcompare -R \
            -o  {params.out_prefix} \
            -r {input.reference} \
            {input.gtfs}"


rule filterUnknownStrand_fasta:
    input:
        "{file}.gtf"
    output:
        "{file}_noUnknownStrand.gtf"
    script:
        "stringtie_filter_fasta.R"


rule sqanti_illumina:
    input:
        "annotation/genome.fa.fai",
        isoforms = "stringtie/illumina/stringtie_illumina_merged_noUnknownStrand.gtf",
        annotation = "annotation/annotation.gtf",
        genome = "annotation/genome.fa"
    output:
        "stringtie/illumina/sqanti/stringtie_illumina_merged_noUnknownStrand_sqanti_report.pdf",
        "stringtie/illumina/sqanti/stringtie_illumina_merged_noUnknownStrand_classification.txt"
    params:
        out_dir = "stringtie/illumina/sqanti"
    shell:
        "sqanti_qc2 -g \
            -d {params.out_dir} \
            {input.isoforms} {input.annotation} {input.genome}"


rule sqanti_flair:
    input:
        "annotation/genome.fa.fai",
        isoforms = "flair/{dataset}/flair.collapse.{dataset}.isoforms.gtf",
        annotation = "annotation/annotation.gtf",
        genome = "annotation/genome.fa"
    output:
        "flair/{dataset}/sqanti/flair.collapse.{dataset}.isoforms_sqanti_report.pdf",
        "flair/{dataset}/sqanti/flair.collapse.{dataset}.isoforms_classification.txt"
    params:
        out_dir = "flair/{dataset}/sqanti"
    shell:
        "sqanti_qc2 -g \
            -d {params.out_dir} \
            {input.isoforms} {input.annotation} {input.genome}"
