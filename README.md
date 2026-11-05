# RNAflowR

**An installable R package for RNA-seq quality control, alignment, gene counting and differential expression.**

RNAflowR runs on **Windows and Linux** from the R console. The same R commands work on both systems; only file paths differ. RNAflowR is research software under development; the initial synthetic-data workflow was tested, while larger real-data validation remains to be completed.

## 1. Install R and RNAflowR (once)

Install a current version of [R](https://cran.r-project.org/). Open **R**, not the Windows Command Prompt or Linux Bash shell.

Run these commands **inside the R console** (`>` prompt):

```r
install.packages(c("BiocManager", "remotes", "pheatmap", "ggplot2", "ggrepel"))
BiocManager::install(c("Rsubread", "edgeR", "limma"), ask = FALSE, update = FALSE)
remotes::install_github("MarioCHMGF/RNAflowR", dependencies = FALSE, upgrade = "never")
library(RNAflowR)
```

**Important:** `install_github()` uses the content actually present on GitHub. The repository must first contain the package files `DESCRIPTION`, `NAMESPACE`, `LICENSE` and `R/` from this release. Until those are uploaded, GitHub installation will not work.

### Optional: add human gene names to the final DEG results

```r
BiocManager::install(c("AnnotationDbi", "org.Hs.eg.db"), ask = FALSE, update = FALSE)
```

These packages are **not required** for read alignment, counting, or statistical analysis. Annotation is handled only by `calcdeg()`; **do not install or load a separate RNAflowR annotation function.**

To update to a newer GitHub version in future, rerun `remotes::install_github(...)`.

## 2. Arrange your files

**Keep reference files and analysis data outside the downloaded package / GitHub repository.** Both Windows and Linux can use this arrangement:

```text
YOUR_HOME/
├── reference_genome_2026/
│   ├── genome.fa                # Reference FASTA (required for indexing)
│   ├── genes.gtf                # OPTIONAL, if not using built-in annotation
│   └── rsubread_index.*         # Generated once by Rsubread
└── RNAflowR_test1/
    ├── files1/                  # Paired-end read 1 (R1) FASTQ.gz
    ├── files2/                  # Paired-end read 2 (R2) FASTQ.gz
    ├── qc/
    ├── bam/
    ├── counts/
    └── results/
```

For example, name paired FASTQs `S01_R1.fastq.gz` and `S01_R2.fastq.gz` and place them in `files1/` and `files2/`, respectively. Single-end datasets only need `files1/`.

Set paths **inside R** as appropriate for your machine:

```r
# Example LINUX:
home_dir <- path.expand("~")
# Example WINDOWS instead: home_dir <- "C:/Users/YourName"
reference_dir <- file.path(home_dir, "reference_genome_2026")
analysis_dir  <- file.path(home_dir, "RNAflowR_test1")
dir.create(reference_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
for (d in c("files1", "files2", "qc", "bam", "counts", "results"))
  dir.create(file.path(analysis_dir, d), showWarnings = FALSE)
setwd(analysis_dir)
```

Copy your FASTQ files to the indicated folders and the matching reference genome FASTA to `reference_dir`. **Do not put large sequence data in GitHub.**

## 3. Index the reference genome (one time per reference)

A **reference genome FASTA** is required for building the Rsubread index. An index made by STAR, BWA or HISAT2 is **not interchangeable** with an Rsubread index. You can reuse an existing Rsubread index without rebuilding.

```r
reference_fasta <- file.path(reference_dir, "genome.fa")  # change to your FASTA filename
index_basename <- file.path(reference_dir, "rsubread_index")

# ONLY if you do not already have this Rsubread index:
prepare_reference(
  reference_fasta = reference_fasta,
  index_basename = index_basename
)
```

The index comprises multiple files with the same basename. When aligning, supply **the basename**, not an individual `.array`, `.tab` or `.files` file.

**No GTF is needed for indexing.** A counting annotation source is chosen separately in Step 6.

## 4. Find and pair your FASTQ samples

The example expects files named `S01_R1.fastq.gz` and `S01_R2.fastq.gz`. Adjust the pattern if your filenames differ.

```r
r1 <- sort(list.files("files1", pattern = "_R1\\.fastq\\.gz$", full.names = TRUE))
r2 <- sort(list.files("files2", pattern = "_R2\\.fastq\\.gz$", full.names = TRUE))
sample_ids <- sub("_R1\\.fastq\\.gz$", "", basename(r1))
r2_ids <- sub("_R2\\.fastq\\.gz$", "", basename(r2))
stopifnot(length(r1) > 0, length(r1) == length(r2),
          !anyDuplicated(sample_ids), setequal(sample_ids, r2_ids))
r2 <- r2[match(sample_ids, r2_ids)]
```

## 5. Check read quality and align

```r
qc <- fastq_qc(c(r1, r2), outdir = "qc")
print(qc)

alignments <- lapply(seq_along(sample_ids), function(i) {
  align_rnaseq(read1 = r1[i], read2 = r2[i],
               index_basename = index_basename,
               sample_id = sample_ids[i], outdir = "bam", nthreads = 8)
})
bam_files <- vapply(alignments, function(x) x$bam, character(1))
```

The built-in quality assessment samples FASTQ base-quality values; it does not trim reads or replace comprehensive FastQC/MultiQC reports. Adjust `nthreads` to your computer.

## 6. Quantify genes: choose ONE counting annotation source

`featureCounts()` requires either **a compatible built-in annotation** or **an external GTF/GFF**. An external `genes.gtf` file is **not obligatory**.

**Option A — Original workflow with built-in hg19 annotation** (use only with a compatible hg19/GRCh37 reference):

```r
counts <- quantify_genes(
  bam_files = bam_files,
  annot_inbuilt = "hg19",
  sample_ids = sample_ids,
  is_paired_end = TRUE,
  strand_specific = 0,
  exclude_sex_chromosomes = FALSE,
  outdir = "counts", nthreads = 8
)
```

**Option B — External GTF/GFF** (necessary for assemblies not supported by built-in annotations, e.g. most plant references):

```r
counts <- quantify_genes(
  bam_files = bam_files,
  annotation_file = file.path(reference_dir, "genes.gtf"),
  sample_ids = sample_ids,
  is_paired_end = TRUE,
  strand_specific = 0,
  exclude_sex_chromosomes = FALSE,
  outdir = "counts", nthreads = 8
)
```

Use **one option only**. The annotation must correspond to the same reference assembly as the FASTA; confirm whether the RNA-seq library is stranded. `is_paired_end = FALSE` for single-end data. `exclude_sex_chromosomes = TRUE` is optional and excludes annotated X/Y genes **after counting and before downstream analysis**. Multimapping reads and multi-overlapping reads are not counted (`FALSE` for both options).

## 7. Provide sample metadata

Create `sampleinfo.csv` in the analysis folder. Samples must have the same identifiers as the count-matrix columns and must be ordered accordingly:

```csv
samples,group
S01,Control
S02,Control
S03,Control
S04,Treatment
S05,Treatment
S06,Treatment
```

Inside R:

```r
sampleinfo <- read.csv("sampleinfo.csv", stringsAsFactors = FALSE)
stopifnot(all(c("samples", "group") %in% names(sampleinfo)),
          !anyDuplicated(sampleinfo$samples),
          setequal(sampleinfo$samples, colnames(counts)))
sampleinfo <- sampleinfo[match(colnames(counts), sampleinfo$samples), , drop = FALSE]
stopifnot(identical(sampleinfo$samples, colnames(counts)))
```

## 8. Run differential expression

The **original CPM filtering method is retained unchanged**; TMM normalization and limma–voom modeling follow. The function makes quality-control figures in the current working directory. Use a dataset with enough genes to produce a heatmap of the 100 most variable genes; tiny toy matrices are not suitable.

```r
setwd(file.path(analysis_dir, "results"))
result <- myrnaseq(
  countdata = counts,
  sampleinfo = sampleinfo,
  samplenum = ncol(counts),
  groupnum = min(table(sampleinfo$group))
)

# Only ONE function is responsible for optional human gene annotations:
deg <- calcdeg(
  dgeObjx = result$dgeObj,
  outcome = sampleinfo$group,
  annotation_db = org.Hs.eg.db::org.Hs.eg.db, # optional; human only
  annotation_keytype = "ENTREZID",            # for built-in hg19 counting IDs
  out_csv = "DEG_annotated.csv"
)

png("volcano.png", width = 1400, height = 1000, res = 150)
print(volcano(deg, Q = TRUE))
dev.off()
```

**To skip annotation**, omit `annotation_db` (and the `annotation_keytype`) in the `calcdeg()` call. For an external human GTF producing Ensembl gene IDs, select `annotation_keytype = "ENSEMBL"` if the identifiers correspond to the database's keys. Version-suffixed IDs (e.g. `.18`) may require conversion to unversioned IDs before matching. For plants and other organisms, annotation can be omitted or performed with a compatible organism-specific `AnnotationDbi` database. `SYMBOL` and `GENENAME` are added **only to the results table**; raw counts, CPM filtering, normalization and model fitting are unaffected. The `myrnaseq()` output `result$degs` remains unannotated.

For an adjusted model, supply up to four aligned covariate vectors to `calcdeg()` as `cov1`, `cov2`, `cov3` and `cov4` (sufficient residual degrees of freedom required).

## 9. Find the output

- `qc/`: per-cycle FASTQ Phred quality summaries and figures.
- `bam/`: aligned BAM files.
- `counts/`: raw gene count matrices, counting statistics and optional X/Y exclusions.
- `results/`: MDS, heatmap, diagnostic figures and DEG results.

## Troubleshooting

- **R command gives a Bash error**: Launch R by typing `R` in Linux Bash, or use the R console on Windows.
- **Package not installed**: Run the installation steps again; `library(RNAflowR)` should work after installation.
- **Reference index not found**: Supply the shared basename of the Rsubread index files. Build only when necessary.
- **Annotation choice error**: Provide exactly one of `annot_inbuilt` or `annotation_file` in `quantify_genes()`; this is a **counting annotation**, distinct from the **optional DEG result annotation** provided by `calcdeg()`.
- **Gene identifiers are unannotated**: Check the organism, identifier type (`ENTREZID` versus `ENSEMBL`) and whether the optional database is installed.
- **Legacy CPM filtering returns undefined threshold**: The original filter depends on genes with exactly 10 counts in each sample; inspect the data rather than silently changing the method.

## Current development status

This package is a reorganized distribution of the existing RNAflowR functions. The initial toy tests are not substitutes for full real-data validation. Scientific results must be reviewed before publication.
