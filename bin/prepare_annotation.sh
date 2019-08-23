#!/bin/bash

# Create folder
mkdir $STORAGE/annotation

# Transcripts
wget -q -O \
    - ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M22/gencode.vM22.transcripts.fa.gz \
    | gunzip > $TX/gencode.vM22.transcripts.fa \

wget -q -O \
    - ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M22/GRCm38.primary_assembly.genome.fa.gz \
    | gunzip > $TX/GRCm38.primary_assembly.genome.fa.gz \

wget -q -O \
    - ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M22/gencode.vM22.annotation.gtf.gz \
    | gunzip > $TX/gencode.vM22.annotation.gtf
