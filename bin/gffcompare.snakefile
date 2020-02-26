rule gffcompare:
    input:
        gtfs = ["stringtie/illumina/stringtie_illumina_merged.gtf",
            "flair/teloprime/flair.collapse.isoforms.gtf",
            "flair/cdna/flair.collapse.isoforms.gtf"],
        reference = "annotation/annotation.gtf"
    output:
        "res/comparisons/gffcompare/all_gffcompare.loci",
        "res/comparisons/gffcompare/all_gffcompare.stats",
        "res/comparisons/gffcompare/all_gffcompare.combined.gtf",
        "res/comparisons/gffcompare/all_gffcompare.tracking"
    params:
        out_prefix = "res/comparisons/gffcompare/all_gffcompare"
    shell:
        "gffcompare -R \
            -o  {params.out_prefix} \
            -r {input.reference} \
            {input.gtfs}"


rule sqanti:
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
