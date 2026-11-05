# Rsubread reference indexing and splice-aware RNA-seq alignment.
# One sample per invocation; for paired-end data pass both mates.
prepare_reference <- function(reference_fasta, index_basename = "reference/index/genome", memory = 8000L) {
  if (!requireNamespace("Rsubread", quietly = TRUE)) stop("Install Rsubread")
  if (!file.exists(reference_fasta)) stop("Missing reference FASTA: ", reference_fasta)
  dir.create(dirname(index_basename), recursive = TRUE, showWarnings = FALSE)
  Rsubread::buildindex(basename = index_basename, reference = reference_fasta,
                       memory = memory)
  invisible(index_basename)
}

align_rnaseq <- function(read1, index_basename, read2 = NULL,
                         sample_id = "sample", outdir = "results/bam", nthreads = 4L,
                         sort_bam = FALSE, phred_offset = 33L, ...) {
  if (!requireNamespace("Rsubread", quietly = TRUE)) stop("Install Rsubread")
  if (length(read1) != 1L || !file.exists(read1)) stop("read1 must be one existing file")
  if (!is.null(read2) && (length(read2) != 1L || !file.exists(read2)))
    stop("read2 must be one existing file or NULL")
  if (!grepl("^[A-Za-z0-9_.-]+$", sample_id)) stop("Invalid sample_id")
  if (!length(Sys.glob(paste0(index_basename, ".*"))))
    stop("Rsubread index not found; call prepare_reference()")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  bam <- file.path(outdir, paste0(sample_id, ".bam"))
  args <- list(index = index_basename, readfile1 = read1,
               type = "rna", input_format = "gzFASTQ", output_format = "BAM",
               output_file = bam, nthreads = as.integer(nthreads),
               phredOffset = phred_offset, sortReadsByCoordinates = sort_bam)
  if (!is.null(read2)) args$readfile2 <- read2
  args <- c(args, list(...))
  alignment_summary <- do.call(Rsubread::align, args)
  if (!file.exists(bam)) stop("Alignment completed without expected BAM output")
  list(bam = bam, alignment_summary = alignment_summary)
}
