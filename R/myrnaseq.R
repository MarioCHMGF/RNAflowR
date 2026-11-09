#' A function to perform RNA seq DEG analysis
#'
#' @param countdata obtained raw counts
#' @param sampleinfo file with mandatory samples(names) and group(what to compare) columns
#' @param samplenum number of all samples
#' @param groupnum number of samples in smallest group
#' @param outdir Directory for saving output figures. Defaults to the current working directory.
#' @return deg dataframe your_object$degs and dgeObj your_object$dgeObj
#' @import edgeR
#' @import limma
#' @import pheatmap
#'
#' @export

myrnaseq <- function(countdata, sampleinfo, samplenum, groupnum,
                     outdir = ".") {

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
 
myCPM <- cpm(countdata)

Fn <- list()

for(i in 1:samplenum){
	temp1 <- data.frame(myCPM[,i],countdata[,i])
	names(temp1) <- c("cpm","count")
	temp2 <- temp1[temp1$count == 10,]
	Fn[[i]] <- mean(temp2$cpm)
}
all <- do.call(rbind, Fn)
useFn <- mean(all[,1])
	
thresh <- myCPM > useFn
keep <- rowSums(thresh) >= groupnum
counts.keep <- countdata[keep,]

dgeObj <- DGEList(counts.keep)

png(file = file.path(outdir, "Barplot_of_library_sizes.png"), width=30, height=15, units="cm", res=300)
barplot(dgeObj$samples$lib.size, names=colnames(dgeObj), las=1)
dev.off()

logcounts <- cpm(dgeObj,log=TRUE)
png(file = file.path(outdir, "Boxplot_of_logCPMs.png"), width=30, height=15, units="cm", res=300)
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=1)
abline(h=median(logcounts),col="blue")
dev.off()

sampleinfo$facs <- factor(sampleinfo$group)

png(file=file.path(outdir, "MDS_dim1_dim2.png"))
levels(sampleinfo$facs)
col.cell <- grDevices::rainbow(nlevels(sampleinfo$facs))[sampleinfo$facs]
data.frame(sampleinfo$group,col.cell)
plotMDS(dgeObj,col=col.cell)
legend("topleft", fill=grDevices::rainbow(nlevels(sampleinfo$facs)), legend=levels(sampleinfo$facs))
title("Status")
dev.off()

var_genes <- apply(logcounts, 1, var)
select_var <- names(sort(var_genes, decreasing=TRUE))[1:100]
highly_variable_lcpm <- logcounts[select_var,]

data_subset <- highly_variable_lcpm

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

data_subset_norm <- t(apply(data_subset, 1, cal_z_score))

my_hclust_gene <- hclust(dist(data_subset_norm), method = "complete")

my_gene_colx <- cutree(my_hclust_gene, k = 2)
my_gene_colx <- data.frame(cluster = ifelse(test = my_gene_colx == 1, yes = "cluster1", no = "cluster2"))

my_sample_col <- data.frame(sample = sampleinfo$group)
row.names(my_sample_col) <- colnames(data_subset)

pheatmap(
  data_subset_norm,
  annotation_row = my_gene_colx,
  annotation_col = my_sample_col,
  cutree_rows = 2,
  cutree_cols = 1,
  filename = file.path(outdir, "High_var_genes.heatmap.png"),
  width = 9,
  height = 12
)

dgeObj <- calcNormFactors(dgeObj)

group <- factor(sampleinfo$group)

design <- model.matrix(~ group)
          
y <- voom(dgeObj, design) 
prefit <- lmFit(y, design)
fit <- eBayes(prefit)
	
png(file = file.path(outdir, "Final_model.png"))
plotSA(fit, main="Final model: Mean−variance trend")
dev.off()

degs <- as.data.frame(topTable(fit, coef=2, confint=T, n = Inf))

newList <- list("dgeObj" = dgeObj, "degs" = degs)

return(newList)
}
