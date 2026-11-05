# Volcano plot: preserves original intent and Q/P options, avoids repeated branches.
volcano <- function(myresults, Q = FALSE, top = 5L, xl = NULL, xu = NULL,
                    logfc_cutoff = 1, p_cutoff = 0.05) {
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("ggrepel", quietly = TRUE)) stop("Install ggplot2 and ggrepel")
  pcol <- if (Q) "adj.P.Val" else "P.Value"
  if (!all(c("logFC", pcol) %in% names(myresults))) stop("Missing required results columns")
  d <- myresults
  d$LOG <- -log10(pmax(d[[pcol]], .Machine$double.xmin))
  d$SIG <- ifelse(!is.na(d[[pcol]]) & d[[pcol]] <= p_cutoff &
                    abs(d$logFC) >= logfc_cutoff, "Significant", "Other")
  if (!"SYMBOL" %in% names(d)) d$SYMBOL <- if ("gene_id" %in% names(d)) d$gene_id else rownames(d)
  d$LABEL <- NA_character_
  candidates <- which(d$SIG == "Significant")
  hits <- head(candidates[order(d[[pcol]][candidates])], top)
  d$LABEL[hits] <- as.character(d$SYMBOL[hits])
  g <- ggplot2::ggplot(d, ggplot2::aes(x = logFC, y = LOG, colour = SIG)) +
    ggplot2::geom_point(alpha = 0.35, size = 1.3, na.rm = TRUE) +
    ggplot2::geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed") +
    ggplot2::geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff), linetype = "dashed") +
    ggrepel::geom_text_repel(data = d[!is.na(d$LABEL), , drop = FALSE],
                             ggplot2::aes(label = LABEL), max.overlaps = Inf) +
    ggplot2::labs(x = "log2 fold change", y = if (Q) "-log10(adjusted p-value)" else "-log10(p-value)", colour = NULL) +
    ggplot2::theme_bw()
  if (!is.null(xl) && !is.null(xu)) g <- g + ggplot2::coord_cartesian(xlim = c(xl, xu))
  g
}
