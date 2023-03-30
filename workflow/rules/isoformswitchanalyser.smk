def get_isa_counts(wildcards):
    if wildcards.dataset == "illumina":
        if wildcards.annotation == "ref":
            return expand("data/quantification/salmon/{barcode}/quant.sf",
                barcode=SAMPLE_INFO_ont["illumina"])
        else:
            return expand("data/quantification/illumina_{{annotation}}/{barcode}/quant.sf",
                barcode=SAMPLE_INFO_ont["illumina"])
    elif wildcards.dataset == "cdna":
        if wildcards.annotation == "ref":
            return expand("data/quantification/cdna/merged/cdna_merged_{barcode}_quant.tsv",
                barcode=SAMPLE_INFO_ont["cdna"])
        else:
            return expand("data/quantification/cdna_{{annotation}}/{barcode}_quant.tsv",
                barcode=SAMPLE_INFO_ont["cdna"])


def get_isa_gtf(wildcards):
    if wildcards.annotation == "ref":
        return "data/annotation/annotation.gtf"
    elif wildcards.annotation == "cdna_flair":
        return "data/reannotation/flair/annotation/cdna_flair.isoforms.gtf"
    elif wildcards.annotation == "illumina_stringtie":
        return "data/reannotation/stringtie/illumina_stringtie_noUnknownStrand.gtf"
    elif wildcards.annotation == "teloprime_stringtie":
        return "data/reannotation/stringtie/teloprime_stringtie_noUnknownStrand.gtf"


def get_isa_transcripts(wildcards):
    if wildcards.annotation == "ref":
        return "data/annotation/transcripts.fa"
    elif wildcards.annotation == "cdna_flair":
        return "data/reannotation/flair/annotation/cdna_flair.isoforms_cleanHeaders.fa"
    elif wildcards.annotation == "illumina_stringtie":
        return "data/reannotation/stringtie/illumina_stringtie_noUnknownStrand.fa"
    elif wildcards.annotation == "teloprime_stringtie":
        return "data/reannotation/stringtie/teloprime_stringtie_noUnknownStrand.fa"


rule isoformswitchanalyser_importData:
    input:
        counts = get_isa_counts,
        sample_info = config["SAMPLE_INFO"],
        gtf = get_isa_gtf,
        transcripts = get_isa_transcripts
    output:
        "data/drimseq/{dataset}_{annotation}/{dataset}_{annotation}_sal.rds",
        "data/drimseq/{dataset}_{annotation}/{dataset}_{annotation}_isoform_AA.fasta",
        res = "data/drimseq/{dataset}_{annotation}/{dataset}_{annotation}_dtu_res.csv.gz"
    params:
        aa_path = "data/drimseq/{dataset}_{annotation}/",
        aa_prefix = "{dataset}_{annotation}_isoform",
        dtu_cutoff = config["dtu_cutoff"], # alpha
        dIF_cutoff = config["dIF_cutoff"],  # min difference
        min_feature_expr = 5,
        min_feature_prop = 0.1,
        min_gene_expr = 10
    wildcard_constraints:
        dataset = "cdna|illumina",
        annotation = "ref|cdna_flair|illumina_stringtie|teloprime_stringtie"
    conda:
        "../envs/r_4.1.2.yaml"
    threads:
        8
    script:
        "../scripts/isoformswitchanalyser_importData.R"


# report
rule isoformswitchanalyser_report:
    input:
        sample_info = config["SAMPLE_INFO"],
        switchList = "data/drimseq/{dataset}_{annotation}/{dataset}_{annotation}_sal.rds",
        pfam = "data/drimseq/{dataset}_{annotation}/{dataset}_{annotation}_isoform_AA.pfam",
        annotation = "data/annotation/annotation_txdb.sqlite",
        telo_tmap = "data/reannotation/stringtie/gffcmp.teloprime_stringtie.gtf.tmap",
        illu_tmap = "data/reannotation/stringtie/gffcmp.illumina_stringtie.gtf.tmap"
    output:
        "res/isoformswitchanalyser/{dataset}_{annotation}.html"
    params:
        dtu_cutoff = config["dtu_cutoff"], # alpha
        dIF_cutoff = config["dIF_cutoff"]  # min difference
    wildcard_constraints:
        dataset = "cdna|illumina"
    conda:
        "../envs/r_4.1.2.yaml"
    script:
        "../scripts/isoformswitchanalyser_report.Rmd"


rule isoformswitchanalyser_all:
    input:
        expand("res/isoformswitchanalyser/{dataset}_{annotation}.html",
            dataset=["cdna", "illumina"],
            annotation=["cdna_flair", "illumina_stringtie", "ref", "teloprime_stringtie"])
