#!/bin/bash

################################################################################
#
# author: Christoph Kiefer
# email: christophak@bmb.sdu.dk
#
# adapter and quality trimming
#
################################################################################

source ../config.sh

mkdir -p $FASTQ/trimmed
mkdir -p ../qc/fastqc/trimmed
mkdir -p ../qc/cutadapt_report

# TruSeq adapters
parallel -j 4 --verbose --eta --results ../qc/cutadapt_report \
    cutadapt \
        -j 6 \
        -q 28 \
        -m 30 \
        -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        -o $FASTQ/trimmed/{}_R1_001_trimmed.fastq.gz \
        -p $FASTQ/trimmed/{}_R2_001_trimmed.fastq.gz \
        $FASTQ/raw/{}_R1_001.fastq.gz \
        $FASTQ/raw/{}_R2_001.fastq.gz \
        :::: ../sample_info/illumina_sample_id.txt

# quality control
ls $FASTQ/trimmed/ | parallel -j 8 --verbose --eta \
    fastqc $FASTQ/trimmed/{} --noextract -o ../qc/fastqc/trimmed

multiqc -f -c .multiqc.conf ../qc -o ../qc

