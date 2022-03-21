# salmon align_reannotation
BEDTOOLS = config["BEDTOOLS"]
MASHMAP = config["MASHMAP"]


rule salmon_prepareDecoys_reannotatedStringtie:
    input:
        genome = "data/annotation/genome.fa",
        annotation = "data/reannotation/stringtie/{dataset}_stringtie_noUnknownStrand.gtf",
        transcripts = "data/reannotation/stringtie/{dataset}_stringtie_noUnknownStrand.fa"
    output:
        "indices/salmon_{dataset}_stringtie/decoy/gentrome.fa",
        "indices/salmon_{dataset}_stringtie/decoy/decoys.txt"
    params:
        outputDir = "indices/salmon_{dataset}_stringtie/decoy"
    threads: 20
    wildcard_constraints:
        dataset = "illumina|teloprime"
    shell:
        "bin/generateDecoyTranscriptome.sh \
            -j {threads} \
            -g {input.genome} \
            -t {input.transcripts} \
            -a {input.annotation} \
            -b {BEDTOOLS} \
            -m {MASHMAP} \
            -o {params.outputDir}"


rule salmon_prepareDecoys_reannotatedFlair:
    input:
        genome = "data/annotation/genome.fa",
        annotation = "data/reannotation/flair/annotation/{dataset}_flair.isoforms.gtf",
        transcripts = "data/reannotation/flair/annotation/{dataset}_flair.isoforms.fa"
    output:
        temp("indices/salmon_{dataset}_flair/decoy/gentrome.fa"),
        temp("indices/salmon_{dataset}_flair/decoy/decoys.txt")
    params:
        outputDir = "indices/salmon_{dataset}_flair/decoy"
    threads: 20
    shell:
        "bin/generateDecoyTranscriptome.sh \
            -j {threads} \
            -g {input.genome} \
            -t {input.transcripts} \
            -a {input.annotation} \
            -b {BEDTOOLS} \
            -m {MASHMAP} \
            -o {params.outputDir}"


rule salmon_index_reannotated:
    input:
        gentrome = "indices/salmon_{dataset}_{method}/decoy/gentrome.fa",
        decoy = "indices/salmon_{dataset}_{method}/decoy/decoys.txt"
    output:
        "indices/salmon_{dataset}_{method}/duplicate_clusters.tsv",
        "indices/salmon_{dataset}_{method}/hash.bin",
        "indices/salmon_{dataset}_{method}/rsd.bin",
        "indices/salmon_{dataset}_{method}/sa.bin",
        "indices/salmon_{dataset}_{method}/txpInfo.bin"
    params:
        outputDir = "indices/salmon_{dataset}_{method}"
    wildcard_constraints:
        method = "flair|stringtie",
        dataset = "illumina|cdna|teloprime"
    threads: 20
    shell:
        "salmon index \
            --gencode \
            -t {input.gentrome} \
            -i {params.outputDir} \
            -d {input.decoy} \
            -p {threads}"


rule salmon_align_reannotated:
    threads: 5
    input:
        "indices/salmon_{dataset}_{method}/duplicate_clusters.tsv",
        "indices/salmon_{dataset}_{method}/hash.bin",
        "indices/salmon_{dataset}_{method}/rsd.bin",
        "indices/salmon_{dataset}_{method}/sa.bin",
        "indices/salmon_{dataset}_{method}/txpInfo.bin",
        fastq_fw = "data/fastq/illumina/trimmed/{sample}_R1_001_trimmed.fastq.gz",
        fastq_rv = "data/fastq/illumina/trimmed/{sample}_R2_001_trimmed.fastq.gz"
    output:
        "data/quantification/illumina_{dataset}_{method}/{sample}/quant.sf"
    params:
        inputDir = "indices/salmon_{dataset}_{method}",
        outputDir = "data/quantification/illumina_{dataset}_{method}/{sample}"
    wildcard_constraints:
        method = "flair|stringtie"
    shell:
        "salmon quant \
            --gcBias \
            --seqBias \
            --numGibbsSamples 25 \
            -i {params.inputDir} \
            -l A \
            -1 {input.fastq_fw} \
            -2 {input.fastq_rv} \
            -p {threads} \
            --validateMappings \
            -o {params.outputDir}"


# minimap
def get_minimapIndexInput(wildcards):
    if wildcards.method == "flair":
        file = "data/reannotation/flair/annotation/{dataset}_flair.isoforms.fa"
    elif wildcards.method == "stringtie":
        file = "data/reannotation/stringtie/{dataset}_stringtie_noUnknownStrand.fa"
    return file


rule minimap_index_reannotation:
    input:
        get_minimapIndexInput
    output:
        "indices/minimap2_{dataset}_{method}/minimap2.mmi"
    threads: 3
    shell:
        "minimap2 \
        -t {threads} \
        -d {output} \
        {input}"


def get_minimapInput(wildcards):
    if wildcards.dataset in ["teloprime", "cdna"]:
        file_name = "data/fastq/{}/{}/{}_q7.fastq.gz".format(
            wildcards.dataset, wildcards.flowcell, wildcards.barcode)
    elif wildcards.dataset == "rna":
        file_name = "data/fastq/rna/{}_q7.fastq.gz".format(wildcards.barcode)
    return file_name


rule minimap_mapReannotation:
    input:
        index = "indices/minimap2_{dataset}_{method}/minimap2.mmi",
        fastq = "data/fastq/cdna/{flowcell}/{barcode}_q7.fastq.gz"
    output:
        sam = temp("data/bam/cdna_{dataset}_{method}/{flowcell}_{barcode}_sort.bam")
    threads: 20
    wildcard_constraints:
        dataset = "teloprime|cdna|illumina",
        type = "flair|stringtie"
    shell:
        "minimap2 \
            -2 \
            -ax map-ont \
            -t {threads} \
            --secondary=no \
            -uf \
            {input.index} \
            {input.fastq} | \
        samtools sort -l 5 -o {output} -O bam -@ 6"


rule samtools_merge2:
    input:
        "data/bam/cdna_{dataset}_{method}/flowcell1_{barcode}_sort.bam",
        "data/bam/cdna_{dataset}_{method}/flowcell2_{barcode}_sort.bam"
    output:
        temp("data/reannotation/cdna_{dataset}_{method}/{barcode}_sort.bam")
    threads: 6
    wildcard_constraints:
        dataset = "teloprime|cdna|illumina",
        method = "stringtie|flair",
        barcode = "barcode.."
    shell:
        "samtools merge -@ {threads} {output} {input}"


rule quantify_minimap_reannotated:
    input:
        "data/reannotation/cdna_{dataset}_{method}/{barcode}_sort.bam"
    output:
        "data/quantification/cdna_{dataset}_{method}/{barcode}_quant.tsv"
    wildcard_constraints:
        dataset = "teloprime|cdna|illumina",
        type = "flair|stringtie"
    script:
        "quantify_minimap.R"


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
    threads:
        8
    script:
        "isoformswitchanalyser_importData.R"


# PFAM
# apt install hmmer2
rule get_pfama:
    output:
        pfam_hmm = "data/pfam/Pfam-A.hmm",
        pfam_hmm_dat = "data/pfam/Pfam-A.hmm.dat",
        pfam_active_site_dat = "data/pfam/active_site.dat"
    shell:
        "wget -q -O \
            - ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz \
            | gunzip > {output.pfam_hmm} && \
        wget -q -O \
            - ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz \
            | gunzip > {output.pfam_hmm_dat} && \
        wget -q -O \
            - ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/active_site.dat.gz \
            | gunzip > {output.pfam_active_site_dat}"


rule hmmpress:
    input:
        pfam_hmm = "data/pfam/Pfam-A.hmm",
        pfam_hmm_dat = "data/pfam/Pfam-A.hmm.dat",
        pfam_active_site_dat = "data/pfam/active_site.dat"
    output:
        multiext("data/pfam/Pfam-A.hmm", ".h3f", ".h3i", ".h3m", ".h3p")
    shell:
        "hmmpress {input.pfam_hmm}"


rule pfam_scan:
    input:
        db = "data/pfam/Pfam-A.hmm",
        fasta = "data/drimseq/{dataset}_{annotation}/{dataset}_{annotation}_isoform_AA.fasta"
    output:
        pfam = "data/drimseq/{dataset}_{annotation}/{dataset}_{annotation}_isoform_AA.pfam",
        tblout = "data/drimseq/{dataset}_{annotation}/{dataset}_{annotation}_isoform_AA.tsv"
    params:
    wildcard_constraints:
        dataset = "cdna|teloprime|illumina|cdna_flair|teloprime_flair"
    threads:
        4
    shell:
        "hmmscan \
            -o {output.pfam} \
            --pfamtblout {output.tblout} \
            --cpu {threads} \
            {input.db} \
            {input.fasta}"


# report
rule isoformswitchanalyser_report:
    input:
        sample_info = config["SAMPLE_INFO"],
        switchList = "data/drimseq/{dataset}_{annotation}/{dataset}_{annotation}_sal.rds"
    output:
        "res/isoformswitchanalyser/{dataset}_{annotation}.html"
    params:
        dtu_cutoff = config["dtu_cutoff"], # alpha
        dIF_cutoff = config["dIF_cutoff"],  # min difference
    script:
        "isoformswitchanalyser_report.Rmd"
