# Description

# Using the pipeline

This analysis makes use of [snakemake](https://snakemake.readthedocs.io/en/stable/).
Run "snakemake --list" for an overview of the available rules.
Run them like this: "snakemake render_ont_gene --cores 20". Snakemake will run
all rules needed to run before running this rule. Adding "-np" will make a dryrun
showing you all the jobs that need to be run in order to fulfill your command.

## Dependencies

### General

* pandoc
* bedtools
* mashmap
* stringtie
* salmon >= 0.14.0
* multiqc
* deeptools
* STAR
* fastqc
* snakemake &gt 5.5.4
* R
* samtools
* pandas
* HTSeq

* libfontconfig1-dev

### R dependencies

R dependencies are handled by renv.

To install the dependencies, start R from the
root folder of the project. If renv is not installed on your computer, it should
install itself into a private library inside the project.

Type the following to restore the library used when this analysis was done:

```R
renv::restore()
```

If the systemfonts package fails to install, you might have to install libcairo2-dev by
running `sudo apt-get libcairo2-dev` first.

# Results files

## Differential gene and transcript expression

### Reports

The folder structure is res/deseq/method/geneOrTxLevel

In every folder, there is a quality control file (..._qc.html) containing sample depths,
sample distances and a PCA plot.

DESeq was run twice, once using log fold change (LFC) shrinking (apeglm) and LFC
cutoff = 1, the other time withouth shrinking and cutoff = 0 (noShrink).

For both DESeq runs, there is a report file (..._report.html), which contains
a table of significantly differentially expressed genes, an MA plot as well as
the results from GO and Reactome analysis. The full DESeq results table is stored
as ..._results.csv.gz.

### Count matrices

There are multiple countmatrices to be found in the res/deseq/method/geneOrTxLevel folders:

* rld: log2(x + 1) level, between sample normalised, variance stabilised. See regularized log (rlog) transformation
[http://dx.doi.org/10.1186/s13059-014-0550-8]. This is the main table for any plotting and correlation analyses etc.

* ntd: log2(x + 1) level, between sample normalised. For genes, that are not expressed in a condition, the rlog transformed
values will not be 0, which might be confusing sometimes. Therfore I add this table.

* tpm: between sample normalised and within sample (between genes) normalised.
(Only for illumina.)

* cts: raw counts, count output from tximport as imported from salmon estimated
counts: **No normalisation at all!** This one should be used for statistical
analyses in DESeq or EdgeR, which perform their own normalisation steps.

#### Using the count files

To avoid any errors reading the count tables back into R (windows/unix, encoding?),
use the following code:

```R
library("readr")
library("dplyr")
cts <- read_csv("res/genelevel/genelevel_cm_cts.csv.gz",
        locale = locale(encoding = "UTF-8"))
```

#### Using the dds object

In the res/deseq/method/geneOrTxLevel folder, there is also a DESeq data set object (dds), that can be imported into R
using readRDS() to be directly analysed in DESeq2 or coerced into an edgeR DGEList object.

```R
# Read the DESeq data set object from disc
dds <- readRDS("res/dds_gencode.vM22_gene.rds")

# Extract the count table and sample info
cts <- counts(dds)
df <- colData(dds)

# create an edgeR DGEList object
library("edgeR")
dgl <- DGEList(counts = cts, samples = df)
```

## Differential transcript usage (DTU)

### DEXSeq

DEXSeq [https://bioconductor.org/packages/release/bioc/vignettes/DEXSeq/inst/doc/DEXSeq.html]
perform a differential exon usage analysis. It compares the fold change of counts falling
into one exon with the counts falling into the whole gene. Therefore it uses a metatranscript
assembled from all transcripts of a gene fullfilling certain expression thresholds.
As it is on exonlevel, it does not depend on correct transcript annotation, but
on correct exon annotation. (It cannot detect differential expression in non annotated
exons.)

It uses feature counts to do the counting. Unfortunatelly, as the reads from the
ONT libraries are not stranded, there is a problem with overlapping genes, causing
many false positives where one exon of a not regulated gene overlaps with one of a
regulated gene.

In the res/dexseq/method folders there is a report (..._dexseq.html) and a heatmap
showing the distribution of the significant exon in the metatranscript (..._heatmap.html).

### DrimSeq

Drimseq [https://bioconductor.org/packages/release/workflows/vignettes/rnaseqDTU/inst/doc/rnaseqDTU.html]
performs differential transcript usage analysis based on the salmon transcript level output.
Therefore it depends on correct transcript annotation. It was run twice, once with the GENCODE
gtf transcriptome and once with a reannotated transcriptome from FLAIR. DrimSeq calculates
a gene and a transcript pvalue and StageR was used to take both in account when
correcting for multiple testing.

The res/dexseq/method folders contain a report (..._drimSeqStageR.html) and a
table with the full StageR output (..._drimSeqStageR.csv).

## Comparisons

The output of all analysis doing some comparisons between the different methods
are in the  res/comparisons folder.
