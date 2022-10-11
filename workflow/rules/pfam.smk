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
    conda:
        "../envs/pfam_1.6.yaml"
    shell:
        "hmmpress {input.pfam_hmm}"


rule pfam_scan:
    input:
        multiext("data/pfam/Pfam-A.hmm", ".h3f", ".h3i", ".h3m", ".h3p"),
        fasta = "data/drimseq/{dataset}_{annotation}/{dataset}_{annotation}_isoform_AA.fasta"
    output:
        "data/drimseq/{dataset}_{annotation}/{dataset}_{annotation}_isoform_AA.pfam"
    params:
        pfam_a_dir = "data/pfam/"
    wildcard_constraints:
        dataset = "cdna|illumina"
    conda:
        "../envs/pfam_1.6.yaml"
    threads:
        5
    shell:
        "pfam_scan.pl \
            -fasta {input.fasta} \
            -dir {params.pfam_a_dir} \
            -outfile {output} \
            -as \
            -cpu {threads}"
