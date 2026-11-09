# RNAflowR

**RNAflowR is an open-source R package for bulk RNA-seq quality assessment, reference-genome indexing, alignment, gene counting, and differential gene expression analysis.**

**Current package version: 0.1.1**  
**Operating systems:** Windows and Linux  
**Source code:** https://github.com/MarioCHMGF/RNAflowR

RNAflowR brings together established tools in a workflow that can be run from the **R console**, with no Bash scripting required. It uses Rsubread for alignment and gene quantification, edgeR for count processing and TMM normalization, and limma–voom for differential expression analysis. The original RNAflowR CPM filtering procedure is retained.

> **Research software:** RNAflowR is under ongoing evaluation. The workflow has been tested with synthetic data and an initial five-sample real-data analysis. Broader benchmarking, additional usability testing, and biological validation are still in progress. Users remain responsible for selecting an appropriate reference genome, library strandedness, sample groups, and statistical design.

---

## 1. Install R and RNAflowR

### 1.1 Install R

Install R from https://cran.r-project.org/. RStudio is optional; RNAflowR also works in the regular R console.

**All commands below are R commands.** On Linux, start an R session by opening a terminal and typing `R`; then enter these commands at the `>` prompt. On Windows, open R or RStudio. Do not paste R commands into Bash or the Windows command prompt.

### 1.2 Install dependencies and RNAflowR

Run the following **once** in R:

```r
install.packages(c("BiocManager", "remotes", "pheatmap", "ggplot2", "ggrepel"))

BiocManager::install(
  c("Rsubread", "edgeR", "limma"),
  ask = FALSE,
  update = FALSE
)

remotes::install_github(
  "MarioCHMGF/RNAflowR",
  dependencies = FALSE,
  upgrade = "never"
)

library(RNAflowR)
packageVersion("RNAflowR")
```

The installed version should be **0.1.1** or newer. Installation is needed only once per R library, although `library(RNAflowR)` must be run in each new R session.

**Optional human gene-name annotation:** If you want gene symbols and descriptions in the final differential-expression table, install:

```r
BiocManager::install(
  c("AnnotationDbi", "org.Hs.eg.db"),
  ask = FALSE,
  update = FALSE
)
```

These optional annotation packages are **not required** for alignment, gene counting, or differential-expression statistics.

### 1.3 Updating RNAflowR

After a new GitHub release, install the current version with:

```r
remotes::install_github(
  "MarioCHMGF/RNAflowR",
  dependencies = FALSE,
  upgrade = "never",
  force = TRUE
)
```

Restart R and load `library(RNAflowR)` afterward.

---

## 2. Prepare the folders

**Keep the reference genome and sequencing files outside the RNAflowR GitHub repository.** Large FASTA, FASTQ, BAM, and index files should not be uploaded to GitHub.

A suggested folder layout:

```text
YOUR_HOME/
├── reference_genome_2026/
│   ├── hg19.fasta                  # Reference genome FASTA
│   ├── hg19_2026_idx.*            # Rsubread index files (if already built)
│   └── genes.gtf                  # OPTIONAL external gene annotation
└── RNAflowR_test1/
    ├── files1/                    # Read 1 FASTQ files
    ├── files2/                    # Read 2 FASTQ files (paired-end only)
    ├── qc/
    ├── bam/
    ├── counts/
    └── results/
```

Create the folders in R:

```r
reference_dir <- path.expand("~/reference_genome_2026")
analysis_dir <- path.expand("~/RNAflowR_test1")

dir.create(reference_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)

for (folder in c("files1", "files2", "qc", "bam", "counts", "results")) {
  dir.create(file.path(analysis_dir, folder), showWarnings = FALSE)
}

setwd(analysis_dir)
getwd()
```

On Windows, `~` refers to the R user's home directory; you may instead set explicit paths, for example `reference_dir <- "C:/RNAflowR_reference"` and `analysis_dir <- "C:/RNAflowR_test1"`. Prefer `/` in paths written in R on both systems.

Use your file manager to copy the reference FASTA into `reference_dir`, and paired FASTQ files into `files1` (R1) and `files2` (R2).

> **Important distinction:** Paired-end sequencing means each biological sample has two read files (R1 and R2). It does **not** imply a paired statistical comparison. Whether samples are matched or independent is a separate aspect of experimental design.

---

## 3. Index the reference genome with Rsubread (once)

To align RNA-seq reads, RNAflowR requires a **reference genome FASTA** and an **Rsubread index built from that FASTA**. An index produced by STAR, HISAT2, or BWA cannot be substituted for an Rsubread index.

### 3.1 If you already have an Rsubread index

Point RNAflowR to its **index basename** (without the `.00.b.array`, `.files`, etc. suffix):

```r
index_basename <- file.path(reference_dir, "hg19_2026_idx")
```

For example, if you have `hg19_2026_idx.00.b.array`, `hg19_2026_idx.00.b.tab`, and other similarly named files, the index basename is `hg19_2026_idx`. **Do not build the index again.**

### 3.2 If you need to build an index

First place your reference genome FASTA in `reference_dir`. Change the FASTA filename to match your file:

```r
reference_fasta <- file.path(reference_dir, "hg19.fasta")
index_basename <- file.path(reference_dir, "hg19_2026_idx")

prepare_reference(
  reference_fasta = reference_fasta,
  index_basename = index_basename
)
```

Index construction may take time and disk space, especially for large genomes, but needs to be done only once per reference/index configuration.

**A GTF/GFF file is not required for genome indexing.** Gene-count annotation is chosen separately in Step 6.

---

## 4. Provide your paired-end FASTQ files

Put one biological sample's R1 and R2 files in the corresponding folders. For example:

```text
files1/SRR8110206_1.fastq     files2/SRR8110206_2.fastq
files1/SRR8110223_1.fastq     files2/SRR8110223_2.fastq
files1/SRR8110224_1.fastq     files2/SRR8110224_2.fastq
files1/SRR8110208_1.fastq     files2/SRR8110208_2.fastq
files1/SRR8110209_1.fastq     files2/SRR8110209_2.fastq
```

**The example below uses the `_1.fastq` / `_2.fastq` naming pattern** commonly produced by `fasterq-dump --split-files`. Edit the filename pattern if your files use `_R1`/`_R2` or another convention. Compressed files may use `.fastq.gz`; check that your installed alignment function supports your input format.

For this example, specify sample IDs once:

```r
sample_ids <- c(
  "SRR8110206", "SRR8110223", "SRR8110224",
  "SRR8110208", "SRR8110209"
)

files1 <- file.path("files1", paste0(sample_ids, "_1.fastq"))
files2 <- file.path("files2", paste0(sample_ids, "_2.fastq"))

stopifnot(all(file.exists(c(files1, files2))))
```

**For your own analysis:** Replace these example identifiers with your samples. `sample_ids` must list the intended samples in the correct order. For compressed files, use `"_1.fastq.gz"` and `"_2.fastq.gz"` instead. For single-end data, only `files1` is needed and `read2` is omitted in alignment.

The five example samples above are a **demonstration**, not a biological study design or a requirement of RNAflowR. Any appropriate sample count can be used, subject to the statistical model.

---

## 5. Assess sequencing quality

```r
qc <- fastq_qc(
  fastq_files = c(files1, files2),
  outdir = "qc",
  nreads = 10000
)

print(qc)
list.files("qc")
```

The default QC assessment samples **10,000 reads from each FASTQ file**. It reports mean Phred quality, Q20/Q30 percentages, and per-cycle summaries and plots. It does **not** examine every sequencing read by default, trim adapters, or replace a comprehensive FastQC/MultiQC assessment.

---

## 6. Align reads to the reference genome

The existing function aligns **one sample at a time**. The following block processes all sample pairs automatically. You do not need to know how the loop works to use it; edit only the paths and number of threads if needed.

```r
alignment_time <- system.time({
  for (i in seq_along(sample_ids)) {
    align_rnaseq(
      read1 = files1[i],
      read2 = files2[i],
      index_basename = index_basename,
      sample_id = sample_ids[i],
      outdir = "bam",
      nthreads = 8
    )
  }
})

alignment_time
list.files("bam", pattern = "\\.bam$")
```

Use a suitable `nthreads` value for your computer. Alignment creates a BAM file for each sample as well as alignment-summary output. The reference index remains in its own folder.

**Single-end data:** Remove the `read2 = files2[i],` line and supply only the single-end read file. **FASTQ format:** The current alignment wrapper configures Rsubread's FASTQ input internally; if it does not recognize a particular compressed/uncompressed input format, check the function and your Rsubread installation before processing a full dataset.

---

## 7. Quantify gene expression with featureCounts

RNAflowR uses `Rsubread::featureCounts()` and returns a **raw gene count matrix**. It saves the count matrix and assignment summaries automatically into `outdir`.

Prepare BAM filenames (no complex R code is needed):

```r
bam_files <- file.path("bam", paste0(sample_ids, ".bam"))
stopifnot(all(file.exists(bam_files)))
```

### IMPORTANT: Determine library strandedness from your sequencing protocol

**Before counting, check the RNA-seq library-preparation manual or the information provided by your sequencing facility.** Do not guess the strandedness based on the experimental groups or file names.

Choose the `strand_specific` value that matches your library:

| Value | Meaning |
|---|---|
| `0` | Unstranded |
| `1` | Forward-stranded |
| `2` | Reverse-stranded |

The default, `strand_specific = 0`, **does not mean your library is necessarily unstranded**. Incorrect strandedness may substantially change read assignment and downstream expression estimates. If the library protocol is unknown, resolve it before treating the biological results as final.

### Option A — Rsubread built-in annotation (no GTF needed)

For a compatible human **hg19/GRCh37** reference, the original RNAflowR approach can use the built-in RefSeq annotation:

```r
counts <- quantify_genes(
  bam_files = bam_files,
  sample_ids = sample_ids,
  annot_inbuilt = "hg19",
  is_paired_end = TRUE,
  strand_specific = 0,  # CHANGE according to your library protocol
  exclude_sex_chromosomes = FALSE,
  outdir = "counts",
  nthreads = 8
)
```

### Option B — External GTF/GFF annotation (when needed)

Use an external annotation for another organism, reference assembly, or annotation release:

```r
counts <- quantify_genes(
  bam_files = bam_files,
  sample_ids = sample_ids,
  annotation_file = file.path(reference_dir, "genes.gtf"),
  is_paired_end = TRUE,
  strand_specific = 0,  # CHANGE according to your library protocol
  exclude_sex_chromosomes = FALSE,
  outdir = "counts",
  nthreads = 8
)
```

Use **exactly one** of `annot_inbuilt` and `annotation_file`. The counting annotation must match the reference genome assembly and its chromosome names. The external GTF/GFF file is **optional** if an appropriate built-in annotation exists.

For single-end data, set `is_paired_end = FALSE`. When `exclude_sex_chromosomes = TRUE`, genes annotated on the selected sex chromosomes are excluded **after quantification and before CPM filtering**, while the full count matrix remains saved. The counting wrapper explicitly uses `countMultiMappingReads = FALSE` and `allowMultiOverlap = FALSE`.

Inspect the resulting matrix and output files:

```r
dim(counts)
head(counts)
colSums(counts)
list.files("counts")
```

The files written by RNAflowR (e.g., `raw_gene_counts.csv`, `featureCounts_assignment_summary.csv`) are generated by the **`quantify_genes()` wrapper**. They are not automatically generated by `featureCounts()` when called directly. The wrapper currently writes these files whenever counting is run.

---

## 8. Create the sample-information table

RNAflowR needs to know which sample belongs to each biological group. Create a small file named **`sampleinfo.csv`** in the analysis directory. The first column contains sample IDs; the second contains groups.

Example for the five samples used above:

```csv
sample,group
SRR8110206,Fibroid
SRR8110223,Fibroid
SRR8110224,Fibroid
SRR8110208,Normal
SRR8110209,Normal
```

A CSV can be created in Excel, LibreOffice Calc, or a text editor. Use the *actual experimental group names and sample IDs* for your dataset.

Then, in R:

```r
sampleinfo <- read.csv("sampleinfo.csv", stringsAsFactors = FALSE)

sampleinfo <- sampleinfo[
  match(colnames(counts), sampleinfo$sample),
  , drop = FALSE
]

stopifnot(
  !anyNA(sampleinfo$sample),
  !anyDuplicated(sampleinfo$sample),
  identical(sampleinfo$sample, colnames(counts))
)

rownames(sampleinfo) <- sampleinfo$sample
```

### Set the direction of your comparison

For a two-group comparison, **the first factor level is the reference group**. If you want positive log2 fold changes to mean higher expression in fibroid tissue:

```r
sampleinfo$group <- factor(
  sampleinfo$group,
  levels = c("Normal", "Fibroid")
)
```

For another dataset, replace these group names with yours. If the study contains paired biological samples, covariates, batches, or more than two groups, the statistical design requires additional consideration; the simple demonstration below is an **unpaired two-group comparison**.

---

## 9. Run CPM filtering, TMM normalization, and limma–voom

RNAflowR's `myrnaseq()` implements the original CPM-based gene filtering procedure, followed by TMM normalization and limma–voom differential-expression analysis. **The filtering procedure has not been changed in version 0.1.1.**

```r
result <- myrnaseq(
  countdata = counts,
  sampleinfo = sampleinfo,
  samplenum = ncol(counts),
  groupnum = min(table(sampleinfo$group)),
  outdir = "results/differential_expression"
)
```

`groupnum` is **the minimum number of samples required to exceed the original CPM threshold**. For this demonstration, it equals the size of the smaller group (2), **not** the number of experimental groups.

The original algorithm calculates a dataset-specific CPM threshold using genes with exactly 10 reads and retains genes above that threshold in at least `groupnum` samples. It has limitations for very small or atypical datasets, including datasets lacking genes with exactly 10 reads; review the filter's suitability for your study.

`myrnaseq()` returns:

```r
result$dgeObj   # Filtered, TMM-normalized DGEList
result$degs     # Differential-expression results (without gene symbols)
```

Check the analysis:

```r
dim(result$dgeObj)
head(result$degs)
list.files("results/differential_expression")
```

**New in v0.1.1:** `outdir` saves the five diagnostic figures in the specified directory, including `High_var_genes.heatmap.png`. The previous intermittent heatmap export problem was addressed by saving through `pheatmap(filename = ...)`.

> The current heatmap displays the 100 most variable retained genes. Very small datasets with too few retained genes or zero-variance genes may require additional handling.

---

## 10. Optional: add gene symbols and gene descriptions

RNAflowR performs **human gene-name annotation only within `calcdeg()`**. There is no second annotation function. This step is optional, and gene symbols are added **only to the final differential-expression results**—not to the raw counts, CPM filter, normalization factors, or model inputs.

For the built-in human hg19 annotation, gene IDs are Entrez identifiers:

```r
degs_annotated <- calcdeg(
  dgeObjx = result$dgeObj,
  outcome = sampleinfo$group,
  annotation_db = org.Hs.eg.db::org.Hs.eg.db,
  annotation_keytype = "ENTREZID",
  out_csv = "results/DEG_annotated.csv"
)

head(degs_annotated)
```

If gene-name annotation is not needed, simply omit `annotation_db` and `annotation_keytype` from the call. For external GTFs with Ensembl IDs, `annotation_keytype = "ENSEMBL"` may be appropriate; version suffixes in Ensembl IDs can require preprocessing. Other organisms need a matching organism-specific annotation database.

**Important implementation detail:** `calcdeg()` currently **refits the limma–voom model** using `result$dgeObj` and the specified `outcome`. For simple two-group comparisons with the same design, results should agree with `myrnaseq()`. For more complex designs, ensure the model specification matches your intended comparison; do not assume this function merely appends labels to an existing results table.

---

## 11. Create a volcano plot

The `volcano()` function returns a ggplot object. The example below uses the adjusted p-value (`Q = TRUE`):

```r
p <- volcano(degs_annotated, Q = TRUE)

print(p)

ggplot2::ggsave(
  filename = "results/volcano.png",
  plot = p,
  width = 8,
  height = 6,
  dpi = 300
)
```

For nominal p-values, use `Q = FALSE`. The plot's default significance and log2 fold-change cutoffs are plotting choices, not automatic proof of biological significance.

---

## 12. Save your analysis for later use

Do not rerun lengthy alignments unnecessarily. Save the main R objects:

```r
saveRDS(
  list(
    sampleinfo = sampleinfo,
    counts = counts,
    result = result,
    degs_annotated = degs_annotated,
    qc = qc
  ),
  "results/RNAflowR_analysis.rds"
)
```

To reopen the analysis in another R session:

```r
saved <- readRDS("results/RNAflowR_analysis.rds")
head(saved$result$degs)
```

The save example assumes you ran optional gene annotation. If you skipped that step, remove `degs_annotated = degs_annotated,` from the list.

For scientific reporting, record the reference assembly and its source, annotation source, FASTQ accessions, library layout and strandedness from the preparation protocol, sample metadata, exact package versions, and computational parameters. You can save package versions using:

```r
writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")
```

---

## 13. Where are the output files?

| Folder | Examples of outputs |
|---|---|
| `qc/` | Per-file QC summary CSV, per-cycle tables and quality plots |
| `bam/` | BAM alignments and Rsubread alignment summaries |
| `counts/` | Raw count matrices, full featureCounts output, assignment statistics, optional chromosome-filter reports |
| `results/differential_expression/` | Library-size barplot, logCPM boxplot, MDS plot, 100-gene heatmap, limma model diagnostic plot |
| `results/` | Annotated DEG CSV, volcano plot, saved analysis, `sessionInfo.txt` |

The exact folders can be changed using each function's `outdir` or file-output arguments.

---

## 14. Troubleshooting

**`bash: syntax error near unexpected token ...`**  
You entered R code in Bash. Start R first (type `R` on Linux) or use R/RStudio on Windows.

**`there is no package called 'RNAflowR'`**  
Complete Step 1 and run `library(RNAflowR)` again.

**Reference index not found**  
Check the `index_basename`. Supply the common index prefix, not one of the individual `.array`, `.tab`, or `.files` paths. An Rsubread index is required.

**FASTQ file not found or wrong pairs**  
Check `sample_ids`, naming conventions, `files1`/`files2`, and `file.exists()`. Paired R1 and R2 files must belong to the same biological sample.

**Low gene-assignment rate**  
Inspect `counts/featureCounts_assignment_summary.csv`. Verify the genome assembly, annotation compatibility, library type, and **strandedness from the library-preparation protocol** before drawing biological conclusions.

**External GTF error**  
Use **either** `annot_inbuilt` **or** `annotation_file`, not both. The reference assembly and chromosome names must agree.

**No annotation symbols**  
Use a suitable organism database and matching identifier key type. `org.Hs.eg.db` is for human annotation only.

**Heatmap missing**  
Check that the installed RNAflowR version is **0.1.1 or newer** and inspect the directory specified by `outdir`. Very small gene matrices may not be suitable for the 100-gene heatmap.

**Unexpected DEG sign (positive/negative logFC)**  
Check the reference factor level in `sampleinfo$group` before running `myrnaseq()`.

---

## 15. Software availability, citation, and development status

RNAflowR is available at **https://github.com/MarioCHMGF/RNAflowR**. It builds on the R/Bioconductor ecosystem and established methods including Rsubread, featureCounts, edgeR, TMM, and limma–voom. When publishing results produced using RNAflowR, also cite the original research publications associated with the underlying methods and the source of any public RNA-seq dataset.

Selected foundational references:

- Liao, Y., Smyth, G. K., & Shi, W. (2013). The Subread aligner: Fast, accurate and scalable read mapping by seed-and-vote. *Nucleic Acids Research*, **41**, e108.
- Liao, Y., Smyth, G. K., & Shi, W. (2014). featureCounts: An efficient general purpose program for assigning sequence reads to genomic features. *Bioinformatics*, **30**, 923–930.
- Liao, Y., Smyth, G. K., & Shi, W. (2019). The R package Rsubread is easier, faster, cheaper and better for alignment and quantification of RNA sequencing reads. *Nucleic Acids Research*, **47**, e47.
- Robinson, M. D., McCarthy, D. J., & Smyth, G. K. (2010). edgeR: A Bioconductor package for differential expression analysis of digital gene expression data. *Bioinformatics*, **26**, 139–140.
- Robinson, M. D., & Oshlack, A. (2010). A scaling normalization method for differential expression analysis of RNA-seq data. *Genome Biology*, **11**, R25.
- Law, C. W., Chen, Y., Shi, W., & Smyth, G. K. (2014). voom: Precision weights unlock linear model analysis tools for RNA-seq read counts. *Genome Biology*, **15**, R29.
- Ritchie, M. E., et al. (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Research*, **43**, e47.
- Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: A practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society: Series B*, **57**, 289–300.

### Version 0.1.1 changes

- Added `outdir` to `myrnaseq()` to keep generated figures together.
- Changed the highly variable gene heatmap to direct PNG export with `pheatmap(filename = ...)`.
- Improved user instructions for library strandedness, reference indexing, file layout, and gene annotation.
- Retained the original CPM filtering, TMM normalization, and limma–voom approach.

**Current scope:** This release supports an accessible, modular RNA-seq processing workflow. Ongoing work includes broader real-data validation, benchmarking, improved error handling, and further simplification of the user interface.
