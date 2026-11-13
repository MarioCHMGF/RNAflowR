# RNAflowR

**A copy–paste-friendly R package for bulk RNA-seq quality checks, alignment, gene counting and differential expression analysis.**

**Version tested:** RNAflowR 0.1.1 · **Platforms:** Linux and Windows · **Code:** https://github.com/MarioCHMGF/RNAflowR

RNAflowR combines Rsubread, featureCounts, edgeR and limma–voom while retaining its author's original CPM expression-filtering procedure. The workflow has been executed with **two six-sample, paired-end public datasets**, one human and one mouse. These are software demonstrations, **not** independent validation of biological conclusions or performance benchmarks.

> **How to use this README:** Copy → paste → run → verify → next step. Commands labeled **R** belong in R/RStudio. Commands labeled **Terminal** belong in a Linux terminal, outside R. Change file locations for your own computer. Research users must confirm reference assembly, annotation, library strandedness and experimental design.

## 1. Install RNAflowR (R; once per R installation)

```r
install.packages(c("BiocManager", "remotes", "pheatmap", "ggplot2", "ggrepel"))
BiocManager::install(c("Rsubread", "edgeR", "limma", "AnnotationDbi"),
                     ask = FALSE, update = FALSE)
remotes::install_github("MarioCHMGF/RNAflowR",
                        dependencies = FALSE, upgrade = "never")
library(RNAflowR)
packageVersion("RNAflowR")
```

Gene-name annotation is optional. For the **human** example, install `org.Hs.eg.db`; for the **mouse** example, install `org.Mm.eg.db`:

```r
BiocManager::install(c("org.Hs.eg.db", "org.Mm.eg.db"),
                     ask = FALSE, update = FALSE)
```

**Update an existing installation:** `remotes::install_github("MarioCHMGF/RNAflowR", force = TRUE, dependencies = FALSE, upgrade = "never")`, then restart R.

## 2. Choose a validated example

| | Human example | Mouse example |
|---|---|---|
| Study | PRJNA498292 / GSE120854 | PRJNA1199025 / GSE284446 |
| Species | *Homo sapiens* | *Mus musculus* |
| Comparison | Fibroid vs normal | MASLD vs standard diet (SD) |
| Samples | 3 + 3 | 3 + 3 |
| Genome | hg19 / GRCh37 | mm10 / GRCm38 |
| featureCounts built-in annotation | `"hg19"` | `"mm10"` |
| Example analysis folder | `~/RNAflowR_test1` | `~/RNAflowR_test2` |

**Both studies used paired-end sequencing; their example statistical comparison is unpaired.** Two FASTQ mates do not imply matched biological cases and controls. The mouse study's sample assignment below is taken from the SRA run metadata.

### 2A. Set up the human experiment (R)

```r
library(RNAflowR)
analysis_dir <- path.expand("~/RNAflowR_test1")
reference_fasta <- path.expand("~/reference_genome_2026/hg19.fasta")
index_basename <- path.expand("~/reference_genome_2026/hg19_2026_idx")
annotation_name <- "hg19"
case_ids <- c("SRR8110206", "SRR8110223", "SRR8110224")
control_ids <- c("SRR8110208", "SRR8110209", "SRR8110210")
case_name <- "Fibroid"
control_name <- "Normal"
```

### 2B. **Or** set up the mouse experiment (R)

Choose this **instead of** 2A if you are processing the mouse data:

```r
library(RNAflowR)
analysis_dir <- path.expand("~/RNAflowR_test2")
reference_fasta <- path.expand("~/reference_genome_mouse/mm10.fa")
index_basename <- path.expand("~/reference_genome_mouse/mm10_idx")
annotation_name <- "mm10"
case_ids <- c("SRR31741942", "SRR31741943", "SRR31741944")
control_ids <- c("SRR31741945", "SRR31741946", "SRR31741947")
case_name <- "MASLD"
control_name <- "SD"
```

**Run all subsequent R steps in the same session** (or rerun the selected setup block when reopening R).

## 3. Create folders (R)

```r
dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
setwd(analysis_dir)
for (folder in c("files1", "files2", "qc", "bam", "counts", "results")) {
  dir.create(folder, showWarnings = FALSE)
}
getwd()
```

Place `_1.fastq` FASTQ files in `files1/`, and their matching `_2.fastq` mates in `files2/`. The filenames must share the same sample ID. Do **not** place FASTQ/BAM/reference files in the GitHub repository.

### Optional: download the exact six public SRA runs (Linux Terminal)

Install the **NCBI SRA Toolkit** to use `fasterq-dump`. These terminal commands download *uncompressed*, split paired-end FASTQs into the current directory. Copy the `_1.fastq` files to `files1/` and `_2.fastq` to `files2/` afterward. Downloads require substantial disk space; `fasterq-dump` also uses temporary storage.

**Human (Terminal):**

```bash
fasterq-dump SRR8110206 --split-files -e 8 -p
fasterq-dump SRR8110223 --split-files -e 8 -p
fasterq-dump SRR8110224 --split-files -e 8 -p
fasterq-dump SRR8110208 --split-files -e 8 -p
fasterq-dump SRR8110209 --split-files -e 8 -p
fasterq-dump SRR8110210 --split-files -e 8 -p
```

**Mouse (Terminal):**

```bash
fasterq-dump SRR31741942 --split-files -e 8 -p
fasterq-dump SRR31741943 --split-files -e 8 -p
fasterq-dump SRR31741944 --split-files -e 8 -p
fasterq-dump SRR31741945 --split-files -e 8 -p
fasterq-dump SRR31741946 --split-files -e 8 -p
fasterq-dump SRR31741947 --split-files -e 8 -p
```

Alternatively, obtain FASTQs through your sequencing facility or another trusted source and copy them with a graphical file manager.

## 4. Obtain/index a matching reference genome (once)

**Mouse mm10 example (Linux Terminal):**

```bash
mkdir -p ~/reference_genome_mouse
cd ~/reference_genome_mouse
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
gunzip mm10.fa.gz
```

**Human hg19:** Download a complete hg19/GRCh37 FASTA from a trusted genome provider and save it as `~/reference_genome_2026/hg19.fasta`. Do not substitute hg38/GRCh38 while using `annot_inbuilt = "hg19"`.

**R (same for either dataset after Step 2A or 2B):**

```r
stopifnot(file.exists(reference_fasta))
if (length(Sys.glob(paste0(index_basename, ".*"))) == 0) {
  prepare_reference(reference_fasta = reference_fasta,
                    index_basename = index_basename)
}
length(Sys.glob(paste0(index_basename, ".*")))  # Expect > 0
```

Index creation can be lengthy. An Rsubread index is not interchangeable with STAR, HISAT2 or BWA indices. **Do not rebuild an index that already exists.** The FASTA assembly must match the selected featureCounts annotation.

## 5. Locate and verify paired FASTQ files (R)

Do not manually type filenames or depend on the order in a directory. These commands discover and match every sample by its ID:

```r
files1 <- sort(list.files("files1", pattern = "_1\\.fastq$", full.names = TRUE))
files2 <- sort(list.files("files2", pattern = "_2\\.fastq$", full.names = TRUE))
sample_ids <- sub("_1\\.fastq$", "", basename(files1))
mate_ids <- sub("_2\\.fastq$", "", basename(files2))

stopifnot(length(files1) == 6, length(files2) == 6,
          identical(sample_ids, mate_ids),
          setequal(sample_ids, c(case_ids, control_ids)),
          all(file.exists(c(files1, files2))))
print(sample_ids)
```

**Expected:** six sample IDs. This example expects **uncompressed** `.fastq` files produced by `fasterq-dump`. If your files are `.fastq.gz`, the discovery patterns and the installed alignment wrapper's supported input format must be checked before using them.

## 6. FASTQ quality control (R)

```r
qc <- fastq_qc(fastq_files = c(files1, files2),
               outdir = "qc", nreads = 10000, min_mean_q = 20)
print(qc)
stopifnot(nrow(qc) == 12, all(qc$passes_mean_q))
list.files("qc")
```

The built-in QC samples the **first 10,000 reads per file**, not every read. It reports mean Phred quality and Q20/Q30 metrics. This is a quick quality check, not a replacement for adapter/content assessment with more comprehensive tools. A QC pass does not establish a correct experimental design.

## 7. Align paired-end reads (R)

This code processes all six samples **sequentially**, with eight alignment threads per sample. It **skips an existing BAM** so a completed sample is not needlessly realigned after an interruption.

```r
for (i in seq_along(sample_ids)) {
  expected_bam <- file.path("bam", paste0(sample_ids[i], ".bam"))
  if (!file.exists(expected_bam)) {
    align_rnaseq(read1 = files1[i], read2 = files2[i],
                 index_basename = index_basename,
                 sample_id = sample_ids[i],
                 outdir = "bam", nthreads = 8)
  }
}

bam_files <- file.path("bam", paste0(sample_ids, ".bam"))
stopifnot(all(file.exists(bam_files)))
print(bam_files)
```

**Important:** `align_rnaseq()` returns a **list** with a BAM path and alignment summary, not a single character string. Do not pass its return value to `vapply(..., FUN.VALUE = character(1))`. The loop above avoids that error. Existing BAMs are skipped by filename; check that they are complete and from the correct reference before reusing them. Inspect the Rsubread mapping summary for each sample.

## 8. Quantify genes with featureCounts (R)

For both worked examples, featureCounts has a **built-in reference annotation**; no external GTF download is necessary:

```r
counts <- quantify_genes(
  bam_files = bam_files,
  sample_ids = sample_ids,
  annot_inbuilt = annotation_name,  # "hg19" or "mm10"
  is_paired_end = TRUE,
  strand_specific = 0,
  exclude_sex_chromosomes = FALSE,
  outdir = "counts",
  nthreads = 8
)

dim(counts)
colSums(counts)
list.files("counts")
```

> **Scientific choice required — library strandedness.** The `0` in the example reproduces the **parameter used for these demonstrations**; it is not proof that the libraries are unstranded. Check the source study or sequencing-provider protocol before making biological claims. Set `strand_specific = 0` for unstranded, `1` for forward-stranded, or `2` for reverse-stranded libraries. An incorrect value can materially distort assignment and DEGs.

FeatureCounts uses `countMultiMappingReads = FALSE` and `allowMultiOverlap = FALSE` in RNAflowR. With `exclude_sex_chromosomes = TRUE`, sex chromosome genes are removed *after counting and before expression filtering*. The full unfiltered gene matrix is still saved. If your assembly is not supported by the built-in annotation, replace `annot_inbuilt = annotation_name` with `annotation_file = "/path/to/matching/genes.gtf"` and use **only one** annotation source.

**Check gene-assignment statistics (R):**

```r
assignment <- read.csv("counts/featureCounts_assignment_summary.csv",
                       check.names = FALSE, row.names = 1)
print(assignment)
stopifnot("Status" %in% names(assignment), "Assigned" %in% assignment$Status)
sample_columns <- setdiff(names(assignment), "Status")
assigned <- as.numeric(assignment[assignment$Status == "Assigned", sample_columns])
total <- colSums(assignment[, sample_columns, drop = FALSE])
print(round(100 * assigned / total, 1))  # % assigned per sample
```

The exported CSV includes a row-number column and a `Status` column; the code deliberately excludes both from the fragment totals.

## 9. Define groups and run differential expression (R)

Group membership is matched **by sample ID**, not alphabetical position. The **control** is the reference group, so positive log2 fold-change indicates higher expression in fibroid/MASLD samples.

```r
sampleinfo <- data.frame(
  sample = colnames(counts),
  group = factor(ifelse(colnames(counts) %in% case_ids,
                        case_name, control_name),
                 levels = c(control_name, case_name)),
  row.names = colnames(counts)
)
stopifnot(setequal(sampleinfo$sample, c(case_ids, control_ids)),
          identical(rownames(sampleinfo), colnames(counts)),
          all(table(sampleinfo$group) == 3))
print(sampleinfo)

result <- myrnaseq(countdata = counts,
                   sampleinfo = sampleinfo,
                   samplenum = ncol(counts),
                   groupnum = min(table(sampleinfo$group)),
                   outdir = "results/DEG")

print(dim(result$degs))
print(sum(result$degs$adj.P.Val < 0.05, na.rm = TRUE))
list.files("results/DEG")
```

**What RNAflowR does:** Its original adaptive CPM threshold is based on the mean CPM of genes having **exactly 10 raw reads in each sample**, averaged across samples. A gene is retained if its CPM exceeds that threshold in at least `groupnum` samples (normally the smaller experimental group size). The retained data are TMM-normalized, transformed using limma–voom, fitted with limma, and evaluated using empirical Bayes moderation and Benjamini–Hochberg FDR adjustment.

`groupnum` is the **minimum number of samples satisfying the CPM threshold**, **not** the number of groups. The simple built-in differential-expression model shown here addresses an **unpaired two-group** comparison without additional covariates. Do not use this example unchanged for paired studies, batch effects, covariate-adjusted designs or multifactor experiments.

RNAflowR saves five diagnostics in `results/DEG/`: `Barplot_of_library_sizes.png`, `Boxplot_of_logCPMs.png`, `MDS_dim1_dim2.png`, `High_var_genes.heatmap.png`, and `Final_model.png`.

## 10. Optional gene symbols and descriptions (R)

Use the organism-specific package corresponding to your chosen example:

```r
annotation_db <- if (annotation_name == "hg19") {
  org.Hs.eg.db::org.Hs.eg.db
} else {
  org.Mm.eg.db::org.Mm.eg.db
}

degs_annotated <- calcdeg(
  dgeObjx = result$dgeObj,
  outcome = sampleinfo$group,
  annotation_db = annotation_db,
  annotation_keytype = "ENTREZID",
  out_csv = "results/DEG_annotated.csv"
)

head(degs_annotated)
```

Both tested built-in annotation choices yield **Entrez gene identifiers**. With other reference annotations (e.g. Ensembl IDs), choose a compatible annotation database and key type. `calcdeg()` **refits the differential-expression model** using the supplied outcome; it does not merely append symbols. For the two example designs, the statistics can be checked for equivalence:

```r
all.equal(result$degs,
          degs_annotated[rownames(result$degs), colnames(result$degs)],
          check.attributes = FALSE)
```

Expected: `TRUE` for the tested simple two-group models. Investigate any difference instead of assuming equivalence.

## 11. Save a volcano plot and the reproducibility bundle (R)

```r
p <- volcano(degs_annotated, Q = TRUE)
ggplot2::ggsave("results/volcano.png", plot = p,
                width = 8, height = 6, dpi = 300)

saveRDS(list(qc = qc, counts = counts, sampleinfo = sampleinfo,
             result = result, annotated_degs = degs_annotated),
        "results/RNAflowR_complete.rds")
writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")
list.files("results")
```

Optional annotation can be skipped; then omit the `degs_annotated` and volcano lines. To restore the saved analysis later: `saved <- readRDS("results/RNAflowR_complete.rds")`.

## 12. What results should I expect from the two examples?

The following results were **observed in RNAflowR 0.1.1** using the named runs, Rsubread built-in annotations, `strand_specific = 0`, no sex chromosome exclusion, and a three-vs-three unpaired comparison. They are **checks for reproducing this specific analysis**, not universal thresholds for RNA-seq data.

| Measure | Human fibroid vs normal | Mouse MASLD vs SD |
|---|---:|---:|
| Quantified genes | 25,702 | 27,179 |
| Genes retained after CPM filter | 18,067 | 13,201 |
| Genes filtered out | 7,635 | 13,978 |
| Significant DEGs, BH FDR < 0.05 | 428 | 4,712 |
| Upregulated in case group | 188 | 2,378 |
| Downregulated in case group | 240 | 2,334 |
| Diagnostic PNGs | 5/5 | 5/5 |

Gene-assignment proportions differed between the experiments (approximately **43.5–60.0% human** vs **84.5–86.7% mouse**). Mapping success and gene-assignment efficiency are **not the same statistic**. Differences can reflect annotation, library properties and analysis settings; they are not grounds on their own for biological conclusions. Likewise, the higher mouse DEG count does **not** establish better performance or higher biological validity.

## 13. Output files and troubleshooting

| Location | What RNAflowR produces |
|---|---|
| `qc/` | FASTQ quality summaries and plots |
| `bam/` | Per-sample BAMs and alignment information |
| `counts/` | `raw_gene_counts.csv/.rds`, `raw_gene_counts_all_chromosomes.csv/.rds`, `featureCounts_assignment_summary.csv`, `featureCounts_full.rds`, sex-chromosome filter reports |
| `results/DEG/` | Five diagnostic PNG figures |
| `results/` | Annotated DEG table, volcano figure, saved R objects and package versions |

**Missing FASTQ / mismatched mates:** Use `print(files1)`, `print(files2)` and the matching assertion in Step 5. Do not infer group identity from sorted file position.

**Reference index missing:** Confirm `index_basename` is the **common prefix**, not an individual `.array` or `.tab` file, and that the reference FASTA really matches the intended assembly.

**Unexpected alignment-loop error:** `align_rnaseq()` returns a list, not one character. Use the `for` loop in Step 7. Avoid rerunning completed alignments simply because a subsequent R call failed.

**Low gene-assignment rate:** Review the assignment summary, reference/annotation compatibility, and **library strandedness**; sequencing mapping percentage is a separate measurement.

**Wrong DEG direction:** Confirm `levels(sampleinfo$group)`; the second factor level represents the case group in these examples.

**Annotation symbols missing:** Confirm `org.Hs.eg.db` for human or `org.Mm.eg.db` for mouse, and a compatible `annotation_keytype`.

**Large data / RAM / disk use:** Downloading FASTQs, building indices and mapping high-depth libraries can consume considerable storage; processing times depend on hardware and dataset size.

## 14. Scope, citation and availability

RNAflowR 0.1.1 has undergone a real-data end-to-end demonstration with human and mouse paired-end datasets, in addition to earlier synthetic checks. These runs demonstrate technical execution and output generation; they are **not** independent accuracy benchmarking against competing RNA-seq workflows, broad multi-platform testing or a substitute for experimental validation.

Source: **https://github.com/MarioCHMGF/RNAflowR**. Cite the underlying Rsubread/Subread, featureCounts, edgeR, TMM and limma–voom methodology when reporting analyses. Foundational papers include Liao *et al.* (2013, 2014, 2019), Robinson *et al.* (2010), Robinson & Oshlack (2010), Law *et al.* (2014), Ritchie *et al.* (2015), and Benjamini & Hochberg (1995). Cite the original publications and accession numbers for the human and mouse datasets if reproducing these examples.

### Version 0.1.1 notes

- `myrnaseq()` writes five diagnostic plots into its user-specified `outdir`, including the high-variance heatmap.
- The original CPM filtering and downstream TMM/limma–voom procedure remain unchanged.
- Mouse and human built-in annotations (`"mm10"`, `"hg19"`) were successfully tested with the same function interfaces.
- `align_rnaseq()` returns a two-part **list** (`bam`, `alignment_summary`), whereas `quantify_genes()` returns an integer gene-by-sample count matrix.
