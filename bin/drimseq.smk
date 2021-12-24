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


# make dmds datasets
# cDNA
rule drimseq_dmdsFromONT:
    input:
        counts = expand("data/quantification/cdna/merged/cdna_merged_{barcode}_quant.tsv", barcode=SAMPLE_INFO_ont["cdna"]),
        txdb = "data/annotation/annotation_txdb.sqlite",
        sample_info = config["SAMPLE_INFO"]
    output:
        "data/drimseq/cdna_ref_dmds.rds"
    script:
        "drimseq_dmdsFromONT.R"


def get_txdb(wildcards):
    if wildcards.method == "flair":
        file =  "data/reannotation/flair/annotation/cdna_flair.isoforms_txdb.sqlite"
    elif wildcards.method == "stringtie":
        if wildcards.dataset == "illumina":
            file = "data/reannotation/stringtie/illumina_stringtie_noUnknownStrand_txdb.sqlite",
        elif wildcards.dataset == "teloprime":
            file = "data/reannotation/stringtie/teloprime_stringtie_noUnknownStrand_txdb.sqlite"
    return file


rule drimseq_dmdsFromONTreannotated:
    input:
        counts = expand("data/quantification/cdna_{{dataset}}_{{method}}/{barcode}_quant.tsv",
            barcode=SAMPLE_INFO_ont["cdna"]),
        txdb = get_txdb,
        sample_info = config["SAMPLE_INFO"]
    output:
        "data/drimseq/cdna_{dataset}_{method}_dmds.rds"
    wildcard_constraints:
        dataset = "teloprime|cdna|illumina",
        type = "flair|stringtie"
    script:
        "drimseq_dmdsFromONT.R"


# illumina
rule tximport_drimseqIlluminaReannotated:
    input:
        salmon_out =
            expand("data/quantification/illumina_{{dataset}}_{{method}}/{sample}/quant.sf",
                sample=SAMPLES_ont),
        txdb = get_txdb,
        sample_info = config["SAMPLE_INFO"]
    output:
        "data/drimseq/illumina_{dataset}_{method}_dmds.rds"
    wildcard_constraints:
        dataset = "teloprime|cdna|illumina",
        type = "flair|stringtie"
    script:
        "drimseq_dmdsFromSalmon.R"


rule tximport_drimseqIllumina:
    input:
        salmon_out = expand("data/quantification/salmon/{sample}/quant.sf", sample=SAMPLES_ont),
        txdb = "data/annotation/annotation_txdb.sqlite",
        sample_info = config["SAMPLE_INFO"]
    output:
        "data/drimseq/illumina_ref_dmds.rds"
    script:
        "drimseq_dmdsFromSalmon.R"


# dtu
rule drimseq_dtu:
    input:
        sample_info = config["SAMPLE_INFO"],
        dmds = "data/drimseq/{file}_dmds.rds"
    output:
        dmds = "data/drimseq/{file}_dmds_dtu.rds"
    params:
        min_feature_expr = 10,
        min_feature_prop = 0.1,
        min_gene_expr = 10
    threads:
        4
    script:
        "drimseq_dtu.R"


rule isoformswitchanalyser_importData:
    input:
        counts = "res/deseq/{dataset}/txlevel/{dataset}_txlevel_cm_cts.csv.gz",
        sample_info = "sample_info/sampleInfo.csv",
        gtf = "annotation/annotation.gtf",
        transcripts = "annotation/transcripts.fa",
        test = "res/drimseq/{dataset}/{dataset}_drimSeqStageR.csv"
    output:
        "res/drimseq/{dataset}/{dataset}_sal.rds",
        "res/drimseq/{dataset}/{dataset}_isoform_AA.fasta"
    params:
        aa_path = "res/drimseq/{dataset}",
        aa_prefix = "{dataset}_isoform"
    wildcard_constraints:
        dataset = "cdna|teloprime|illumina"
    script:
        "drimseq_isoformswitchanalyser_importData.R"


rule clean_fasta_ids:
    input:
        "{file}.fa"
    output:
        temp("{file}_clean.fa")
    shell:
        "sed 's/_[^_]*$//' {input} > {output}"


rule isoformswitchanalyser_importData_flair:
    input:
        counts = "flair/{dataset}/flair_{dataset}_counts_matrix.tsv",
        sample_info = "sample_info/sampleInfo.csv",
        gtf = "flair/{dataset}/flair.collapse.{dataset}.isoforms.gtf",
        transcripts = "flair/{dataset}/flair.collapse.{dataset}.isoforms_clean.fa",
        test = "res/drimseq/{dataset}_flair/{dataset}_flair_drimSeqStageR.csv"
    output:
        "res/drimseq/{dataset}_flair/{dataset}_flair_sal.rds",
        "res/drimseq/{dataset}_flair/{dataset}_flair_isoform_AA.fasta"
    params:
        aa_path = "res/drimseq/{dataset}_flair",
        aa_prefix = "{dataset}_flair_isoform"
    wildcard_constraints:
        dataset = "cdna|teloprime"
    script:
        "drimseq_isoformswitchanalyser_importData_flair.R"


rule signalP:
    input:
        "res/drimseq/{dataset}/{dataset}_isoform_AA.fasta"
    output:
        "res/drimseq/{dataset}/{dataset}_isoform_AA_summary.signalp5"
    params:
        out_prefix = "res/drimseq/{dataset}/{dataset}_isoform_AA"
    wildcard_constraints:
        dataset = "cdna|teloprime|illumina|cdna_flair|teloprime_flair"
    threads:
        40
    shell:
        "signalp \
            -fasta {input} \
            -format short \
            -org euk \
            -plot none \
            -prefix {params.out_prefix}"


rule targetP:
    input:
        "res/drimseq/{dataset}/{dataset}_isoform_AA.fasta"
    output:
        "res/drimseq/{dataset}/{dataset}_isoform_AA_summary.targetp2"
    params:
        out_prefix = "res/drimseq/{dataset}/{dataset}_isoform_AA"
    wildcard_constraints:
        dataset = "cdna|teloprime|illumina|cdna_flair|teloprime_flair"
    threads:
        40
    shell:
        "targetp \
            -fasta {input} \
            -format short \
            -org non-pl \
            -plot none \
            -prefix {params.out_prefix}"


def get_nt_fasta(wildcards):
    if wildcards.dataset in ["teloprime", "illumina", "cdna"]:
        fasta = "annotation/transcripts.fa"
    elif wildcards.dataset == "teloprime_flair":
        fasta = "flair/teloprime/flair.collapse.teloprime.isoforms.fa"
    elif wildcards.dataset == "cdna_flair":
        fasta = "flair/cdna/flair.collapse.cdna.isoforms.fa"
    return fasta


rule CPAT:
    input:
        fasta = get_nt_fasta,
        hex = "data/cpat/Mouse_Hexamer.tsv",
        logit = "data/cpat/Mouse_logitModel.RData"
    output:
        "res/drimseq/{dataset}/{dataset}_cpat.ORF_prob.best.tsv"
    params:
        out_prefix = "res/drimseq/{dataset}/{dataset}_cpat"
    wildcard_constraints:
        dataset = "cdna|teloprime|illumina|cdna_flair|teloprime_flair"
    threads:
        40
    shell:
        "cpat.py \
            -g {input.fasta} \
            -o {params.out_prefix} \
            -d {input.logit} \
            -x {input.hex} && \
         rm {params.out_prefix}.ORF_seqs.fa {params.out_prefix}.ORF_prob.tsv \
            {params.out_prefix}.no_ORF.txt {params.out_prefix}.r"


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
        multiext("data/pfam/Pfam-A.hmm", ".h3f", ".h3i", ".h3m", ".h3p"),
        fasta = "res/drimseq/{dataset}/{dataset}_isoform_AA.fasta"
    output:
        "res/drimseq/{dataset}/{dataset}_isoform_AA.pfam"
    params:
        pfam_a_dir = "data/pfam/"
    wildcard_constraints:
        dataset = "cdna|teloprime|illumina|cdna_flair|teloprime_flair"
    threads:
        5
    shell:
        "pfam_scan \
            -fasta {input.fasta} \
            -dir {params.pfam_a_dir} \
            -outfile {output} \
            -as \
            -cpu {threads}"


def get_axis(wildcards):
    if wildcards.dataset == "illumina":
        axis = "illumina"
    else:
        axis = "ont"
    return axis


rule drimseq_report:
    input:
        dmds = "res/drimseq/{dataset}/{dataset}_dmds_dtu.rds",
        res = "res/drimseq/{dataset}/{dataset}_drimSeqStageR.csv",
        signalP = "res/drimseq/{dataset}/{dataset}_isoform_AA_summary.signalp5",
        targetP = "res/drimseq/{dataset}/{dataset}_isoform_AA_summary.targetp2",
        pfam = "res/drimseq/{dataset}/{dataset}_isoform_AA.pfam",
        cpat = "res/drimseq/{dataset}/{dataset}_cpat.ORF_prob.best.tsv",
        biomaRt_gene = "annotation/biomaRt_gene.rds",
        biomaRt_tx = "annotation/biomaRt_tx.rds",
        switchList = "res/drimseq/{dataset}/{dataset}_sal.rds"
    output:
        "res/drimseq/{dataset}/{dataset}_drimSeqStageR.html"
    params:
        axis = get_axis
    wildcard_constraints:
        dataset = "teloprime|cdna|illumina|cdna_flair|teloprime_flair"
    script:
        "drimseq_stagerAnalysis.Rmd"
