# Make bigwigs
rule index_BAM:
    input:
        "BAM/{sample}_Aligned.sortedByCoord.out.bam"
    output:
        "BAM/{sample}_Aligned.sortedByCoord.out.bam.bai"
    shell:
        "samtools index {input}"

rule make_bigwigs:
    input:
        bam = "BAM/{sample}_Aligned.sortedByCoord.out.bam",
        bai = "BAM/{sample}_Aligned.sortedByCoord.out.bam.bai"
    output:
        fw = "BW/{sample}_fw.bw",
        rv = "BW/{sample}_rv.bw"
    threads: 10
    shell:
        "bamCoverage \
            -b {input.bam} \
            -o {output.fw} \
            --filterRNAstrand forward \
            -p {threads} \
            --effectiveGenomeSize 2652783500  \
            --normalizeUsing BPM && \
        bamCoverage \
            -b {input.bam} \
            -o {output.rv} \
            --filterRNAstrand reverse \
            -p {threads} \
            --effectiveGenomeSize 2652783500 \
            --normalizeUsing BPM"

rule make_hub:
    input:
        expand("BW/{sample}_fw.bw", sample=SAMPLES),
        expand("BW/{sample}_rv.bw", sample=SAMPLES),
        sample_info = "sample_info/sampleInfo.csv"
    output:
        "nanoporeibat_hub/hub/hub.txt",
        "nanoporeibat_hub/hub/genomes.txt",
        "nanoporeibat_hub/hub/mm10/trackDb.txt",
        expand("nanoporeibat_hub/bw/{sample}_fw.bw", sample=SAMPLES),
        expand("nanoporeibat_hub/bw/{sample}_rv.bw", sample=SAMPLES)
    params:
        url = "http://bioinformatik.sdu.dk/solexa/webshare/christoph/nanoporeibat_hub"
    shell:
        "bin/generateUCSChub.R {input.sample_info} {params.url}"
