# Corrected optional-covariate design; original model remains limma-voom.
# outcome is the factor/continuous predictor of interest; covariates can be named columns
# in sampleinfo or explicit vectors aligned to the count matrix columns.
calcdeg <- function(dgeObjx, outcome, cov1 = NULL, cov2 = NULL, cov3 = NULL,
                    cov4 = NULL, sampleinfo = NULL, annotation_db = NULL,
                    annotation_keytype = "ENTREZID", out_csv = NULL) {
  if (!requireNamespace("limma", quietly = TRUE)) stop("Install limma")
  n <- ncol(dgeObjx$counts)
  resolve <- function(x, label) {
    if (is.null(x)) return(NULL)
    if (is.character(x) && length(x) == 1L && !is.null(sampleinfo) && x %in% names(sampleinfo))
      x <- sampleinfo[[x]]
    if (length(x) != n || anyNA(x)) stop(label, " must be complete and length ", n)
    x
  }
  outcome <- resolve(outcome, "outcome")
  if (is.character(outcome)) outcome <- factor(outcome)
  covs <- list(cov1 = cov1, cov2 = cov2, cov3 = cov3, cov4 = cov4)
  covs <- covs[!vapply(covs, is.null, logical(1))]
  dat <- data.frame(outcome = outcome)
  for (nm in names(covs)) dat[[nm]] <- resolve(covs[[nm]], nm)
  if (is.factor(dat$outcome) && nlevels(dat$outcome) < 2L) stop("outcome must contain at least two levels")
  design <- stats::model.matrix(stats::reformulate(names(dat)), data = dat)
  if (qr(design)$rank < ncol(design)) stop("Design matrix is rank deficient")
  if (nrow(design) <= ncol(design)) stop("Insufficient residual degrees of freedom")
  coefs <- grep("^outcome", colnames(design), value = TRUE)
  if (is.numeric(outcome)) coefs <- "outcome"
  if (length(coefs) != 1L) stop("Choose a two-level outcome; multi-level contrasts need an explicit design")
  y <- limma::voom(dgeObjx, design)
  fit <- limma::eBayes(limma::lmFit(y, design))
  results <- as.data.frame(limma::topTable(fit, coef = coefs, confint = TRUE, n = Inf))
  results$gene_id <- rownames(results)
  results <- results[, c("gene_id", setdiff(names(results), "gene_id")), drop = FALSE]
  if (!is.null(annotation_db)) {
    if (!requireNamespace("AnnotationDbi", quietly = TRUE)) stop("Install AnnotationDbi")
    mapped <- AnnotationDbi::select(annotation_db, keys = results$gene_id,
                      keytype = annotation_keytype, columns = c("SYMBOL", "GENENAME"))
    mapped <- mapped[!duplicated(mapped[[annotation_keytype]]), , drop = FALSE]
    idx <- match(results$gene_id, mapped[[annotation_keytype]])
    results$SYMBOL <- mapped$SYMBOL[idx]
    results$GENENAME <- mapped$GENENAME[idx]
  }
  if (!is.null(out_csv)) {
    dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
    write.csv(results, out_csv, row.names = FALSE)
  }
  results
}
