# library('DESeq2')
library('edgeR')
library('testit')

peform_dge <- function(dge.object, minReplicatesForReplace) 
{
  dge.object <- DESeq2::DESeq(dge.object, test = "Wald",
                              minmu = 0.5, 
                              minReplicatesForReplace = minReplicatesForReplace, 
                              parallel = FALSE)
  
  # how well does the model describe the data?
  # plotDispEsts(dge.object)
  
  result_names <- DESeq2::resultsNames(dge.object)
  
  
  # Shrink lof2fc
  # cite: Zhu, A., Ibrahim, J.G., Love, M.I. (2018) Heavy-tailed prior distributions for 
  # sequence count data: removing the noise and preserving large differences. 
  # Bioinformatics. 10.1093/bioinformatics/bty895 
  res.lfc <- DESeq2::lfcShrink(dge.object, coef=tail(result_names, n=1), type="apeglm", 
                               lfcThreshold = 0, format = c("DataFrame"))
  res.lfc <- res.lfc[order(res.lfc$padj), ]
  #   remove NA's from result
  res.lfc <- na.omit(res.lfc)
  
  # Without log2fc shirnkage 
  dge.result <- DESeq2::results(dge.object, 
                                name = tail(result_names, n=1), # USE LAST TERM OF RESULTS NAME
                                alpha = 0.05, pAdjustMethod = "BH", 
                                tidy = TRUE, format = c("DataFrame"))
  
  return(list(res.lfc, dge.result, dge.object))
}



#' @description Run DGE analysis for general experiments (with multiple factors) 
def.run_general_edger <- function(y.obj, design.matrix)  
{
  # Perform DGE analysis 
  message("Perform DGE analysis for general experiments")
  
  y.obj <- edgeR::calcNormFactors(y.obj)
  
  # Estimate dispersion - robust against outliers
  y.obj <- edgeR::estimateDisp(y.obj, design.matrix, robust=TRUE)
  # edgeR::plotBCV(y.obj)
  
  fit <- edgeR::glmQLFit(y.obj, design.matrix, robust=TRUE) 
  # edgeR::plotQLDisp(fit)
  
  # Differential expression test - takes by default last column in design
  qlf <- edgeR::glmQLFTest(fit)
  
  # Multiple correction test
  fdr.value <- p.adjust(qlf$table$PValue, method="BH")
  qlf$table$padj <- fdr.value
  
  df.genefeatures <- y.obj$genes
  testit::assert(all(row.names(qlf$table) == df.genefeatures$hgnc_symbol))
  dge.result = cbind(
    qlf$table, base::data.frame(
      df.genefeatures[, c('description', 'entrezid', 'ENSEMBL')]))
  
  #   sort summary by adjusted p-val
  message("Create summary Dataframe")
  dge.result <- dge.result[order(dge.result$padj), ]
  
  dge.result <- dge.result %>% dplyr::rename(log2FoldChange = logFC)
  dge.result$hgnc_symbol <- row.names(dge.result)
  
  return(list(dge.result, qlf))
}

