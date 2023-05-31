__author__ = "Christoph Kiefer"
__email__ = "christophak@bmb.sdu.dk"


import pandas as pd


configfile: "conf/config.yaml"


# Sample IDs
SAMPLE_INFO = pd.read_csv(config["SAMPLE_INFO"], sep=',')
SAMPLE_INFO = SAMPLE_INFO.set_index("sample_id")
SAMPLE_INFO_illumina = SAMPLE_INFO[SAMPLE_INFO["illumina"].notnull()]
SAMPLES = SAMPLE_INFO_illumina["illumina"].tolist()
SAMPLE_INFO_ont = SAMPLE_INFO[SAMPLE_INFO["ont"].notnull()]
SAMPLES_ont = SAMPLE_INFO_ont["illumina"].tolist()
BARCODES = SAMPLE_INFO_illumina[SAMPLE_INFO_illumina["ont"].notnull()]["ont"].tolist()


rule all:
    input:
        "multiqc_trimmed.html",
        "res/browserTracks/Adcy3.pdf",
        "res/browserTracks/Cars2.pdf",
        "res/comparisons/comparisons_countsPCA.html",
        "res/comparisons/comparisons_coverage.html",
        "res/comparisons/comparisons_dgeDte.html",
        "res/comparisons/comparisons_dtu.html",
        "res/comparisons/comparisons_feature_detection.html",
        "res/comparisons/comparisons_go.html",
        "res/comparisons/comparisons_quantification_correlation.html",
        "res/comparisons/comparisons_quantification_correlation_normalised.html",
        "res/comparisons/comparisons_quantification_correlationWithinSamples.html",
        "res/comparisons/comparisons_quantification_correlationWithinSamples_normalised.html",
        "res/comparisons/comparisons_readLengths_bam_genome.html",
        "res/comparisons/comparisons_readLengths_bam_transcriptome.html",
        "res/comparisons/comparisons_readLengths_fastq.html",
        "res/comparisons/comparisons_reannotation.html",
        "res/qpcr/qpcr_dtu_validation.html",
        "res/qpcr/qpcr_temperature_effect.html"


# Include other rules
include: "bin/annotation.smk"

include: "bin/process_illumina.smk"
include: "bin/process_nanopore.smk"

include: "bin/rseqc.smk"

include: "bin/comparisons.smk"

include: "bin/deseq.smk"

include: "workflow/rules/stringtie.smk"
include: "workflow/rules/flair.smk"
include: "bin/reannotation_quant.smk"

include: "workflow/rules/sqanti2.smk"
include: "workflow/rules/gffcompare.smk"
include: "bin/browserTracks.smk"

include: "workflow/rules/isoformswitchanalyser.smk"

include: "workflow/rules/qpcr.smk"
include: "workflow/rules/pfam.smk"

include: "workflow/rules/compare_quantification.smk"
include: "workflow/rules/compare_featureDetection.smk"
include: "workflow/rules/compare_reannotationDTU.smk"
include: "workflow/rules/compare_dge.smk"

