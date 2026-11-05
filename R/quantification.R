# Gene-level quantification using Rsubread::featureCounts.
# Optionally exclude sex-chromosome genes from the gene count matrix BEFORE
# downstream CPM filtering, normalization, and differential expression.
# No modification is made to the existing CPM filtering algorithm.

#' Remove genes annotated on sex chromosomes from a featureCounts matrix.
#'
#' @param counts Numeric gene-by-sample count matrix, with gene IDs as row names.
#' @param annotation The `annotation` data.frame returned by featureCounts,
#'   containing `GeneID` and `Chr` columns.
#' @param sex_chromosomes Reference-specific chromosome identifiers. Default
#'   recognizes human-style X/Y and chrX/chrY. Override for other organisms.
#' @return A list with `counts`, `removed_gene_ids`, and `removed_n`.
exclude_sex_chromosome_genes <- function(
    counts, annotation,
    sex_chromosomes = c("X", "Y", "chrX", "chrY")) {
  if (!is.matrix(counts)) stop("counts must be a gene-by-sample matrix")
  if (is.null(rownames(counts)) || anyDuplicated(rownames(counts)))
    stop("counts must have unique gene ID row names")
  if (!is.data.frame(annotation) ||
      !all(c("GeneID", "Chr") %in% colnames(annotation)))
    stop("featureCounts annotation must include GeneID and Chr columns")
  if (!is.character(sex_chromosomes) || !length(sex_chromosomes) ||
      anyNA(sex_chromosomes) || any(!nzchar(sex_chromosomes)))
    stop("sex_chromosomes must contain chromosome names")
  if (anyDuplicated(as.character(annotation$GeneID)))
    stop("Annotation has duplicated GeneID values")

  # featureCounts sometimes lists multiple chromosomes for a gene using ';'.
  # Use exact token matching to avoid inadvertently matching e.g. chrX_random.
  chromosome_tokens <- strsplit(as.character(annotation$Chr), "[;,]", perl = TRUE)
  located_on_sex_chr <- vapply(chromosome_tokens, function(tokens) {
    any(trimws(tokens) %in% sex_chromosomes)
  }, logical(1))
  chromosome_by_id <- setNames(located_on_sex_chr,
                               as.character(annotation$GeneID))
  if (!all(rownames(counts) %in% names(chromosome_by_id)))
    stop("Some count-matrix gene IDs are missing from featureCounts annotation")
  remove <- unname(chromosome_by_id[rownames(counts)])
  list(counts = counts[!remove, , drop = FALSE],
       removed_gene_ids = rownames(counts)[remove],
       removed_n = sum(remove))
}

#' Count gene-level reads and optionally exclude sex-chromosome genes.
#'
#' @param annotation_file Optional external GTF/GFF filename. Supply exactly one of
#'   `annotation_file` and `annot_inbuilt`.
#' @param annot_inbuilt Optional Rsubread built-in annotation, e.g. 'hg19'.
#'   The original author's script used `annot.inbuilt='hg19'`.
#' @param exclude_sex_chromosomes Logical; FALSE by default. When TRUE, genes
#'   mapped to any of `sex_chromosomes` are excluded from returned counts.
#' @param sex_chromosomes Names matching the reference annotation's Chr column.
#'   For example c("X", "Y") or c("chrX", "chrY").
#' Multimapping reads and multi-overlapping reads are not counted.
#' @return Count matrix for downstream analysis (filtered if requested).
quantify_genes <- function(bam_files, annotation_file = NULL, sample_ids = NULL,
                           is_paired_end = FALSE, strand_specific = 0L,
                           outdir = "results/counts", nthreads = 4L,
                           feature_type = "exon", gene_attribute = "gene_id",
                           annot_inbuilt = NULL,
                           exclude_sex_chromosomes = FALSE,
                           sex_chromosomes = c("X", "Y", "chrX", "chrY")) {
  if (!requireNamespace("Rsubread", quietly = TRUE)) stop("Install Rsubread")
  if (!length(bam_files) || any(!file.exists(bam_files))) stop("Missing BAM files")
  has_file <- !is.null(annotation_file) && length(annotation_file) == 1L &&
    !is.na(annotation_file) && nzchar(annotation_file)
  has_inbuilt <- !is.null(annot_inbuilt) && length(annot_inbuilt) == 1L &&
    !is.na(annot_inbuilt) && nzchar(annot_inbuilt)
  if (has_file == has_inbuilt)
    stop("Supply exactly one of annotation_file (GTF/GFF) or annot_inbuilt (e.g. 'hg19')")
  if (has_file && !file.exists(annotation_file)) stop("Missing GTF/GFF annotation")
  if (is.null(sample_ids)) sample_ids <- sub("\\.bam$", "", basename(bam_files), ignore.case = TRUE)
  if (length(sample_ids) != length(bam_files) || anyDuplicated(sample_ids))
    stop("Provide unique sample IDs matching BAM files")
  if (!strand_specific %in% 0:2) stop("strand_specific must be 0, 1, or 2")
  if (!is.logical(exclude_sex_chromosomes) || length(exclude_sex_chromosomes) != 1L ||
      is.na(exclude_sex_chromosomes))
    stop("exclude_sex_chromosomes must be TRUE or FALSE")

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  fc_args <- list(files = bam_files, useMetaFeatures = TRUE,
                  countMultiMappingReads = FALSE, allowMultiOverlap = FALSE,
                  isPairedEnd = is_paired_end, countReadPairs = is_paired_end,
                  strandSpecific = as.integer(strand_specific),
                  nthreads = as.integer(nthreads))
  if (has_inbuilt) {
    fc_args$annot.inbuilt <- annot_inbuilt
  } else {
    fc_args$annot.ext <- annotation_file
    fc_args$isGTFAnnotationFile <- TRUE
    fc_args$GTF.featureType <- feature_type
    fc_args$GTF.attrType <- gene_attribute
  }
  fc <- do.call(Rsubread::featureCounts, fc_args)
  counts <- fc$counts
  colnames(counts) <- sample_ids
  storage.mode(counts) <- "integer"

  # Always keep the complete featureCounts output so filtering is auditable.
  saveRDS(fc, file.path(outdir, "featureCounts_full.rds"))
  saveRDS(counts, file.path(outdir, "raw_gene_counts_all_chromosomes.rds"))
  write.csv(data.frame(gene_id = rownames(counts), counts, check.names = FALSE),
            file.path(outdir, "raw_gene_counts_all_chromosomes.csv"), row.names = FALSE)
  write.csv(fc$stat, file.path(outdir, "featureCounts_assignment_summary.csv"))

  removed <- character(0)
  if (exclude_sex_chromosomes) {
    result <- exclude_sex_chromosome_genes(counts, fc$annotation,
                                          sex_chromosomes = sex_chromosomes)
    counts <- result$counts
    removed <- result$removed_gene_ids
  }
  write.csv(data.frame(gene_id = removed),
            file.path(outdir, "excluded_sex_chromosome_genes.csv"), row.names = FALSE)
  write.csv(data.frame(
      exclude_sex_chromosomes = exclude_sex_chromosomes,
      chromosomes = if (exclude_sex_chromosomes) paste(sex_chromosomes, collapse = ",") else "",
      total_genes = nrow(fc$counts), excluded_genes = length(removed),
      remaining_genes = nrow(counts)),
      file.path(outdir, "chromosome_filter_summary.csv"), row.names = FALSE)
  saveRDS(counts, file.path(outdir, "raw_gene_counts.rds"))
  write.csv(data.frame(gene_id = rownames(counts), counts, check.names = FALSE),
            file.path(outdir, "raw_gene_counts.csv"), row.names = FALSE)
  counts
}
