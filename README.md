# RNAflowR — development starter

A modular extension of the author's earlier `myrnaseq`, `calcdeg`, and `volcano` R functions. **Development code: not yet validated on a complete FASTQ dataset.** Retain the original repository and its attribution/history.

## Existing core

- `R/myrnaseq.R`: original CPM filtering **unchanged**, downstream exploratory analysis and limma-voom; MDS filename and a few plotting issues corrected.
- `R/calcdeg.R`: corrected adjustable model with up to four covariates; supports a species-specific `AnnotationDb` object.
- `R/volcano.R`: simplified plotting function; both nominal and adjusted p-values.

## New upstream modules

- `R/quality_control.R`: `fastq_qc()` extracts a sample of FASTQ Phred scores via `Rsubread::qualityScores`, writes per-cycle metrics and plots. **This is limited to sampled base quality and is NOT full FastQC/MultiQC, trimming, contamination, adapter detection, or a guarantee of sequencing suitability.** Use FastQC/MultiQC for comprehensive QC.
- `R/alignment.R`: `prepare_reference()` builds a Rsubread index; `align_rnaseq()` generates BAM files from single/paired-end FASTQ.
- `R/quantification.R`: `quantify_genes()` runs gene-level `featureCounts` with either Rsubread built-in hg19 annotation or a compatible GTF/GFF, outputs counts, full object and assignment metrics. Optional `exclude_sex_chromosomes=TRUE` excludes annotated X/Y genes **after counting and before CPM filtering**, preserving the full unfiltered counts for audit.

## Prerequisites

R with Bioconductor `Rsubread`, `edgeR`, `limma`, `AnnotationDbi` (optional species annotation), `pheatmap`, `dendextend`, plus `ggplot2`, `ggrepel` and original `myrnaseq` dependencies. Install relevant packages from CRAN/Bioconductor.

## Quick start

```r
source("R/quality_control.R")
source("R/alignment.R")
source("R/quantification.R")
source("R/myrnaseq.R")
source("R/calcdeg.R")
source("R/volcano.R")

# One reference build, then align each sample. FASTQ paths are examples.
prepare_reference("reference/genome.fa", "reference/index/genome")
fastq_qc(c("fastq/S1_R1.fastq.gz", "fastq/S1_R2.fastq.gz"))
a1 <- align_rnaseq("fastq/S1_R1.fastq.gz", "reference/index/genome",
                   read2="fastq/S1_R2.fastq.gz", sample_id="S1")
# Repeat for S2, S3, ... and pass all BAMs at once for quantification.
# counts <- quantify_genes(c(a1$bam, a2$bam, ...), "reference/genes.gtf",
#                          sample_ids=c("S1", "S2", ...), is_paired_end=TRUE,
#                          exclude_sex_chromosomes=TRUE,
#                          sex_chromosomes=c("X", "Y", "chrX", "chrY"))
# sampleinfo <- data.frame(samples=colnames(counts), group=c(...))
# result <- myrnaseq(counts, sampleinfo, ncol(counts), min(table(sampleinfo$group)))
# deg <- calcdeg(result$dgeObj, outcome=sampleinfo$group,
#                out_csv="results/degs.csv")
# print(volcano(deg, Q=TRUE))
```

**Important:** Confirm the reference genome, GTF chromosome naming, paired-end setting, strandedness and sample order before running. For human data, the annotation database can be `org.Hs.eg.db::org.Hs.eg.db` if row IDs are Entrez IDs; Ensembl IDs require a matching `annotation_keytype`. No species annotation is required to run the statistical modeling.

## Known methodological boundary

The exact pre-existing CPM logic uses genes with **exactly 10 reads** in each sample to estimate a threshold, as instructed by the project owner. This may produce undefined values when such genes are absent; the method has deliberately **not** been altered. Validate against an expected dataset before publication.

`myrnaseq()` currently performs an unadjusted two-group `limma-voom` fit and `calcdeg()` performs a separate optional covariate-adjusted fit. `myrnaseq()` expects sample columns in exactly the same order as `sampleinfo` rows. `plotMDS()` produces MDS, not PCA.

This repository starter has been prepared without an R interpreter in the authoring environment; full runtime validation and biological benchmarking remain required.

## Optional sex-chromosome exclusion

`quantify_genes(..., exclude_sex_chromosomes=FALSE)` **keeps all chromosomes** by default. Set it to `TRUE` to remove genes assigned by featureCounts to any chromosome in `sex_chromosomes` (defaults: `X`, `Y`, `chrX`, `chrY`). Use exact chromosome names from the GTF/GFF; these names may differ in nonhuman species and references. Gene entries with more than one annotated chromosome are excluded if **any** of those chromosomes is in the list. Mitochondrial chromosomes are **not** excluded.

All raw counts are saved as `raw_gene_counts_all_chromosomes.*` and the counts for downstream analysis as `raw_gene_counts.*`. The removal list (`excluded_sex_chromosome_genes.csv`) and summary (`chromosome_filter_summary.csv`) make filtering reproducible. This happens **before** `myrnaseq()`, whose pre-existing CPM filtering algorithm has not been modified. If you run `myrnaseq()` directly with an externally supplied count matrix, chromosome filtering does **not** occur automatically; use the helper `exclude_sex_chromosome_genes(counts, featureCounts_annotation, sex_chromosomes=...)` first.

## Original built-in hg19 annotation compatibility

The original manual script used `featureCounts(bam.files, annot.inbuilt="hg19", isPairedEnd=TRUE)`.
The modular workflow now supports this by specifying `annot_inbuilt="hg19"` **instead** of `annotation_file`. Do not supply both.

```r
counts <- quantify_genes(
  bam_files = bam.files,
  annot_inbuilt = "hg19",
  is_paired_end = TRUE,
  exclude_sex_chromosomes = TRUE
)
```

The function returns a matrix with **gene IDs as row names** and **one column per BAM sample**, so no sample-number-dependent column selection is required. The complete matrix, filtered matrix, and exclusions are saved. The original hard-coded `m[, c(2:9)]` is intentionally not reproduced.

For reproducible new analyses and nonhuman organisms, a versioned reference assembly and matching external annotation are preferable.


### Counting policy
`quantify_genes()` explicitly sets `countMultiMappingReads = FALSE` and `allowMultiOverlap = FALSE`, so multimapping alignments and multi-overlapping reads do not contribute to gene counts. Sex-chromosome exclusion remains optional. Original CPM threshold calculations are unchanged.
