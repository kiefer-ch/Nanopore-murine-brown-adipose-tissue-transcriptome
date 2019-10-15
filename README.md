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
* samtools
* pandas
* HTSeq

* libcairo2-dev
* libfontconfig1-dev

R dependencies are handled by packrat. Run "snakemake packrat_init" to restore
the libraries needed for this analysis. This will create a private package
library for this project, ensuring the same package versions are used as were
used originally.

## Results files

### Differential gene and transcript expression

I have done the differential expression analyses twice: once with only the samples,
that have been sequenced on the ont machine too (_ont), and once with all samples that
have been sequenced on the illumina machine (_all).

There are multiple countmatrices to be found in the res/genelevel_... and
res/txlevel_... folders:

* rld: log2 level, between sample normalised, variance stabilised. See regularized log (rlog) transformation
[http://dx.doi.org/10.1186/s13059-014-0550-8]. This is the main table for any plotting and correlation analyses etc.

* ntd: log2 level, between sample normalised. For genes, that are not expressed in a condition, the rlog transformed
values will not be 0, which might be confusing sometimes. Therfore I add this table.

* tpm: between sample normalised and within sample (between genes) normalised.

* cts: raw counts, count output from tximport as imported from salmon estimated
counts: **No normalisation at all!** This one should be used for statistical
analyses in DESeq or EdgeR, which perform their own normalisation steps.

In the data/ folder, there is also the DESeq data set object (dds), that can be imported into R
using readRDS() to be directly analysed.
