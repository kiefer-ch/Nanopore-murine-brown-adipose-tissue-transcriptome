# nanoPoreIBAT

## Description

## Using the pipeline

This analysis makes use of [snakemake](https://snakemake.readthedocs.io/en/stable/).
Run "snakemake --list" for an overview of the available rules. 
Run them like this: "snakemake render_ont_gene --cores 20". Snakemake will run 
all rules needed to run before running this rule. Adding "-np" will make a dryrun 
showing you all the jobs that need to be run in order to fulfill your command.

## Dependencies

* pandoc
* bedtools
* mashmap
* stringtie
* salmon >= 0.14.0
* multiqc
* deeptools
* STAR
* fastqc
* snakemake
* R
            
R dependencies are handled by packrat. Run "snakemake packrat_init" to restore 
the libraries needed for this analysis. This will create a private package 
library for this project, ensuring the same package versions are used as were
used originally.
