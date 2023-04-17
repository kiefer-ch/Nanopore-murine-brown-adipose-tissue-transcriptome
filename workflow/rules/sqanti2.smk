rule filterUnknownStrand_fasta:
    input:
        "{file}.gtf"
    output:
        "{file}_noUnknownStrand.gtf"
    conda:
        "../envs/r_4.1.2.yaml"
    script:
        "../scripts/gffcompare_filterFasta.R"


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
