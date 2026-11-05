# Basic FASTQ base-quality QC with Rsubread. Does not trim reads or replace FastQC/MultiQC.
# Sampled reads are used to avoid loading full FASTQ files into memory.
fastq_qc <- function(fastq_files, outdir = "results/qc", nreads = 10000L,
                     phred_offset = 33L, min_mean_q = 20) {
  if (!requireNamespace("Rsubread", quietly = TRUE)) stop("Install Bioconductor package Rsubread")
  if (!is.character(fastq_files) || !length(fastq_files) || any(!file.exists(fastq_files)))
    stop("fastq_files must contain existing FASTQ paths")
  if (length(nreads) != 1L || is.na(nreads) || nreads < 1L) stop("nreads must be positive")
  if (!phred_offset %in% c(33L, 64L)) stop("phred_offset must be 33 or 64")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  summary_rows <- vector("list", length(fastq_files))
  for (i in seq_along(fastq_files)) {
    file <- fastq_files[[i]]
    q <- Rsubread::qualityScores(filename = file, offset = phred_offset, nreads = as.integer(nreads))
    if (!is.matrix(q) || !nrow(q)) stop("No quality scores recovered from: ", file)
    q <- matrix(as.numeric(q), nrow = nrow(q), ncol = ncol(q))
    q[!is.finite(q) | q < 0] <- NA_real_
    if (all(is.na(q))) stop("No usable Phred scores recovered from ", file)
    per_cycle <- data.frame(cycle = seq_len(ncol(q)),
                            mean_phred = colMeans(q, na.rm = TRUE),
                            pct_q20 = colMeans(q >= 20, na.rm = TRUE) * 100,
                            pct_q30 = colMeans(q >= 30, na.rm = TRUE) * 100)
    sample_label <- sprintf("sample_%02d", i)
    write.csv(per_cycle, file.path(outdir, paste0(sample_label, "_per_cycle.csv")), row.names = FALSE)
    png(file.path(outdir, paste0(sample_label, "_quality.png")), width = 1600, height = 900, res = 160)
    plot(per_cycle$cycle, per_cycle$mean_phred, type = "l", lwd = 2,
         xlab = "Read position (bp)", ylab = "Mean Phred quality",
         main = basename(file), ylim = c(0, max(40, per_cycle$mean_phred, na.rm = TRUE)))
    abline(h = c(20, 30), lty = 2, col = "grey55")
    dev.off()
    summary_rows[[i]] <- data.frame(file = normalizePath(file), reads_sampled = nrow(q),
      mean_phred = mean(q, na.rm = TRUE), pct_q20 = mean(q >= 20, na.rm = TRUE) * 100,
      pct_q30 = mean(q >= 30, na.rm = TRUE) * 100,
      passes_mean_q = mean(q, na.rm = TRUE) >= min_mean_q)
  }
  report <- do.call(rbind, summary_rows)
  write.csv(report, file.path(outdir, "fastq_quality_summary.csv"), row.names = FALSE)
  report
}
