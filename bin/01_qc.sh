#!/bin/bash

################################################################################
#
# author: Christoph Kiefer
# email: christophak@bmb.sdu.dk
#
# inital quality control
#
################################################################################

# config.sh contains the path to the folder holding the fastq files
source ../config.sh

mkdir -p ../qc/fastqc/raw

ls $FASTQ/raw/ | parallel -j 4 --verbose --eta \
    fastqc $FASTQ/raw/{} --noextract -o ../qc/fastqc/raw

multiqc -f ../qc/fastqc/raw -o ../qc 

