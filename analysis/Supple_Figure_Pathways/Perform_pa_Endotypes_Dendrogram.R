# Script to perform DGE analysis between branches of a node
rm(list = ls(all=TRUE))

library('xlsx')
library(openxlsx)
library(data.table)
library("dendextend")
# Pathway analysis libraries
library(ReactomePA)
library('fgsea')
library('DOSE')
library('org.Hs.eg.db')

general_path = '/Volumes/CH__data/Projects/Eyerich_AG_projects/BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer/analysis'
source(file.path(general_path, 'Molecular_subtypes/r_scripts/Pathway_analysis/ORA_GSEA_plots.R'))


#' @description
#' Extract significantly differentially expressed genes (DEGs) from the DGE result 
#' @return returns a list containing a dataframe with DEGs, EntrezIDs of 
#' the DEGs, and their HGNC symbol  
#'  
get_significantgenes <- function(df.dge_results, padj.cut, lfc_factor, op) 
{
  # II. Get significantly differentially expressed genes
  df.sig <- df.dge_results[df.dge_results$padj < padj.cut & 
                             !is.na(df.dge_results$padj) & 
                             op(df.dge_results$log2FoldChange, lfc_factor), ]
  degenes.sig <- df.sig$entrezid
  degenes.sig <- as.character(na.exclude(degenes.sig))       
  
  degenes.sig.hgnc_symbol <- df.sig$hgnc_symbol
  degenes.sig.hgnc_symbol <- as.character(na.exclude(degenes.sig.hgnc_symbol))  
  
  return(list(df.sig, degenes.sig, degenes.sig.hgnc_symbol))
}

#' @description
#' Rank genes based on their Log2FC 
#' @return list of including two lists of ranked genes one time with HGNC symbol 
#' and the other with EntrezIDs 
#' 
do_rank_genes <- function(df.dge_results) 
{
  # II. Rank all genes based on their fold change
  ranked_genes.genesymbol <- df.dge_results$log2FoldChange
  names(ranked_genes.genesymbol) <- df.dge_results$hgnc_symbol
  ranked_genes.genesymbol <- sort(ranked_genes.genesymbol, decreasing = T)
  
  # with entrezID as names
  ranked_genes.entrezid <- df.dge_results$log2FoldChange
  names(ranked_genes.entrezid) <- df.dge_results$entrezid
  ranked_genes.entrezid <- sort(ranked_genes.entrezid, decreasing = T)
  
  return(list(ranked_genes.genesymbol, ranked_genes.entrezid))
}


# 0. Parameters 
obs.name <- "healthysamp_clusterlabels"
opt.method = 'cosine'
opt.linkage_method = 'ward.D2'

selected.database <- 'ReactomePA'
fdr.value <- 0.05
padj.cut_degs <- 0.05
l2fc.factor <- 1
minGSSize <- 10

multitest.method <- 'BH'
pval.cut <- 0.05

date.data <- '2023-08-29'
general.dir <- file.path(general_path, 'output')
input.dir <- file.path(general.dir, 'DGE_analysis', 'Endotypes_Dendrogram')
output.dir <- file.path(general.dir, 'Pathways', 'Endotypes_Dendrogram', Sys.Date())
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

# 1. Load data
df.grouped.counts <- read.csv(file.path(general_path, 'output/Dendrogram/2023-09-13/Normed_Dataframe_mean_gene_selection.csv'))
sorting <- order(as.integer(gsub("E", "\\1", df.grouped.counts$cluster)))
# Reorder rows by row name
df.grouped.counts <- df.grouped.counts[ sorting, ]
cluster.names <- df.grouped.counts$cluster
row.names(df.grouped.counts) <- seq(1, dim(df.grouped.counts)[1])
df.grouped.counts$cluster <- NULL

# 2. Create dendrogram
Matrix <- as.matrix(df.grouped.counts)
res.dist <- Matrix / sqrt(rowSums(Matrix * Matrix))
res.dist <- res.dist %*% t(res.dist)
res.dist <- as.dist(1 - res.dist)
Z <- hclust(d=res.dist, method=opt.linkage_method)
tree <- ape::as.phylo(Z)

# 3. Create Pathway analysis pairs
df.pairs <- copy(Z$merge)
for (ind in seq(dim(Z$merge)[1]))
{
  sum.rows <- sum(Z$merge[ind, ] < 0)
  if (sum.rows <= 1)
  {
    # replace column(s) with row entry
    pos.values <- Z$merge[ind, ]
    # Which column contains the positive value?
    ind.pos.value <- which(Z$merge[ind, ] > 0)
    
    if (length(ind.pos.value) > 1) 
    {
      df.pairs[ind, 1] <- paste(df.pairs[pos.values[1], ], collapse = ", ")
      df.pairs[ind, 2] <- paste(df.pairs[pos.values[2], ], collapse = ", ")
    } else 
    {
      df.pairs[ind, ind.pos.value] <- paste(df.pairs[pos.values[ind.pos.value], ], collapse = ", ")
    }
  } 
}


# 4. Run Pathway enrichment one vs one
message('Pathway enrichment analysis of L skin biopsies within branches of Dendrogram')
hash.dict <- hash::hash()
# Save upregulated gene lists to list
DEGs_entrezID <- vector(mode = "list", length = dim(df.pairs)[1])
DEGs_hgnc <- vector(mode = "list", length = dim(df.pairs)[1])
DEGs.ranked_genes <- vector(mode = "list", length = dim(df.pairs)[1])
ora_list <- vector(mode = "list", length = dim(df.pairs)[1])
gsea_list <- vector(mode = "list", length = dim(df.pairs)[1])
obj.gsea_list <- vector(mode = "list", length = dim(df.pairs)[1])


for (row_number in seq(1, dim(df.pairs)[1])) 
{
  df.comb.temp <- df.pairs[as.integer(row_number), ]
  val.1 <- as.integer(unlist(strsplit(df.comb.temp[1],", "))) * -1
  val.2 <- as.integer(unlist(strsplit(df.comb.temp[2],", "))) * -1
  
  val.1 <- paste("E", val.1, sep='')
  val.2 <- paste("E", val.2, sep='')
  
  
  message(paste0(paste(val.1, collapse = '_'), ' vs ',
                 paste(val.2, collapse = '_')))
  
  tryCatch(
    {
      design.name = 'Age_Sex_batchID_subtypes'
      comparison.name <- paste(paste(val.1, collapse = '_'), ' vs ',
                               paste(val.2, collapse = '_'), design.name, sep = "__")

      
      # 4. Load data
      df.dge <- read.csv(file.path(input.dir, date.data, 
                                   paste0('DEGs_', comparison.name, '.csv')))
      
      
    }, error=function(e)
    {
      design.name = 'Age_Sex_subtypes'
      comparison.name <- paste(paste0(paste(val.1, collapse = '_'), ' vs ',
                                      paste(val.2, collapse = '_')), 
                               design.name, sep = "__")
      
      message('Switch to design without correcting for batchID')
      # 4. Load data
      df.dge <- read.csv(file.path(input.dir, date.data, 
                                   paste0('DEGs_', comparison.name, '.csv')))
    })
  
  # Save parameters used in this study
  df.parameters <- data.frame(
    Parameters=c("DGE comparison", "", 'Log2FC_cut_DEGs', 'padj-value_cut_DEGs', 'FDR_cut_PA',
                 'p-value_cut_PA', 'Multitesting_method', 'minGSSize', 'database', '',
                 'Preselected PAs', 'File_preselected_PAs', 'sheetName'),
    Values=c(file.path(input.dir, date.data, 
                       paste0('DEGs_', comparison.name, '.csv')), "", l2fc.factor, padj.cut_degs, 
             fdr.value, pval.cut, multitest.method, minGSSize, 'Reactome',
             '', 'no', '', 'CH_Reactome'))
  
  # ========================= Prepare Data ========================= #
  # I. a) Add signed p.adj values
  df.dge$sign.padj <- -sign(df.dge$log2FoldChange) * log10(df.dge$padj)
  
  # I. b) Extract upregulated genes in lesion Endotype skin
  df.up_degenes = get_significantgenes(df.dge_results = df.dge, padj.cut = padj.cut_degs,
                                        lfc_factor = l2fc.factor, op = `>`) 
  df.up = df.up_degenes[[1]]
  # Save significant genes with entrezID
  DEGs_entrezID[[paste(val.1, collapse = '_')]] <- df.up_degenes[[2]]
  # Save significant genes with HGNC Symbol
  DEGs_hgnc[[paste(val.1, collapse = '_')]] <- df.up_degenes[[3]]
  
  # I. c) Extract down regulated genes in lesion Endotype skin
  df.down_degenes = get_significantgenes(df.dge_results = df.dge, padj.cut = padj.cut_degs,
                                       lfc_factor = -l2fc.factor, op = `<`) 
  df.down = df.down_degenes[[1]]
  # Save significant genes with entrezID
  DEGs_entrezID[[paste(val.2, collapse = '_')]] <- df.down_degenes[[2]]
  # Save significant genes with HGNC Symbol
  DEGs_hgnc[[paste(val.2, collapse = '_')]] <- df.down_degenes[[3]]
  
  # I. d) Rank genes based on their fold change
  ranked_genes.gs_enID.up = do_rank_genes(df.dge_results = df.up)
  ranked_genes.gs_enID.down = do_rank_genes(df.dge_results = df.down)
  
  # I. e) Rank all genes based on their fold change
  ranked_genes.gs_enID.all = do_rank_genes(df.dge_results = df.dge)
  ranked_genes.enID.all <- ranked_genes.gs_enID.all[[2]]
  ranked_genes.gs.all <- ranked_genes.gs_enID.all[[1]]
  DEGs.ranked_genes[[comparison.name]] <- ranked_genes.gs.all
  
  # ========================= Perform Pathway analysis ========================= #
  # Perform ORA of upregulated DEx using ReactomePA DB
  pa.L_up <- ReactomePA::enrichPathway(
    gene=names(ranked_genes.gs_enID.up[[2]]), pvalueCutoff=pval.cut, readable=T, universe = df.dge$entrezid)
  ora_list[[paste(val.1, collapse = '_')]] <- pa.L_up@result[
    pa.L_up@result$p.adjust < fdr.value, ]$Description
  pa.L_down <- ReactomePA::enrichPathway(
    gene=names(ranked_genes.gs_enID.down[[2]]), pvalueCutoff=pval.cut, 
    readable=T, universe = df.dge$entrezid)
  ora_list[[paste(val.2, collapse = '_')]] <- pa.L_down@result[
    pa.L_down@result$p.adjust < fdr.value, ]$Description
  
  # GSEA
  gsea.all <- ReactomePA::gsePathway(
    ranked_genes.enID.all, organism = "human", eps = 0,
    pvalueCutoff=pval.cut, pAdjustMethod=multitest.method, verbose=FALSE,
    minGSSize = minGSSize, by = "fgsea", seed=TRUE)
  gsea.all <- DOSE::setReadable(gsea.all, 'org.Hs.eg.db', 'ENTREZID')
  gsea_list[[comparison.name]] <- gsea.all@result[gsea.all@result$p.adjust < pval.cut, ]$Description
  obj.gsea_list[[comparison.name]] <- gsea.all@result
  
  
  message("Save to xlsx file")
  xlsx::write.xlsx(pa.L_up@result, file.path(
    output.dir,  paste('ORA__', paste(val.1, collapse = '_'), 
                       "__comparison_", comparison.name, '.xlsx', 
                       sep='')), sheetName="result")
  xlsx::write.xlsx(df.parameters, file = file.path(
    output.dir,  paste('ORA__', paste(val.1, collapse = '_'), 
                       "__comparison_", comparison.name, '.xlsx', 
                       sep='')), sheetName="info", append=TRUE)
  xlsx::write.xlsx(pa.L_down@result, file.path(
    output.dir,  paste('ORA__', paste(val.2, collapse = '_'), 
                       "__comparison_", comparison.name, '.xlsx',
                       sep='')), sheetName="result")
  xlsx::write.xlsx(df.parameters, file = file.path(
    output.dir,  paste('ORA__', paste(val.2, collapse = '_'), 
                       "__comparison_", comparison.name, '.xlsx', 
                       sep='')), sheetName="info", append=TRUE)
  
  xlsx::write.xlsx(gsea.all@result, file.path(
    output.dir, paste('GSEA__', comparison.name, '.xlsx', sep='')), sheetName="result")
  xlsx::write.xlsx(df.parameters, file = file.path(
    output.dir,  paste('GSEA__', comparison.name, '.xlsx', sep='')), sheetName="info", append=TRUE)
  
  
  # 5. Plots
  # ORA
  if (!is.null(nrow(pa.L_up)))
  {
    if (nrow(pa.L_up) > 1 & length(unique(pa.L_up@result$ID)) > 1)
    {
      def.plot.enricher_barplot(enricher.result=pa.L_up, showCategories=6,
                                method='ORA: Reactome',
                                title=paste(selected.database, paste(val.1, collapse = '_'),
                                            "__comparison_", comparison.name,
                                            'ORA_barplot.pdf', sep = "_"),
                                width=5, height=4, output.dir=output.dir, func_colours=NULL)
      
      # def.pathways.cnetplot(enrich.result=pa.L_up,
      #                       method='ORA: Reactome',
      #                       entrezid_log2fc=ranked_genes.gs.all, showCategories=3,
      #                       title=paste(selected.database,  paste(val.1, collapse = '_'),
      #                                   'ORA_cnetplot.pdf', sep = "_"),
      #                       width=12, height=6, output.dir=output.dir)
    }
  }
  
  if (!is.null(nrow(pa.L_down)))
  {
    if (nrow(pa.L_down) > 1 & length(unique(pa.L_down@result$ID)) > 1)
    {
      def.plot.enricher_barplot(enricher.result=pa.L_down, showCategories=6,
                                method='ORA: Reactome',
                                title=paste(selected.database, paste(val.2, collapse = '_'),
                                            "__comparison_", comparison.name,
                                            'ORA_barplot.pdf', sep = "_"),
                                width=5, height=4, output.dir=output.dir, func_colours=NULL)
      
      # def.pathways.cnetplot(enrich.result=pa.L_down,
      #                       method='ORA: Reactome',
      #                       entrezid_log2fc=ranked_genes.gs.all, showCategories=3,
      #                       title=paste(selected.database,  paste(val.2, collapse = '_'),
      #                                   'ORA_cnetplot.pdf', sep = "_"),
      #                       width=12, height=6, output.dir=output.dir)
    }
  }
  
  # GSEA
  if (!is.null(nrow(gsea.all)))
  {
    if (nrow(gsea.all) > 1 & length(unique(gsea.all@result$ID)) > 1)
    {
      # def.pathways.cnetplot(enrich.result=gsea.all,
      #                       method='GSEA: Reactome',
      #                       entrezid_log2fc=ranked_genes.gs.all, showCategories=3,
      #                       title=paste(selected.database, comparison.name,
      #                                   'GSEA_cnetplot.pdf', sep = "_"),
      #                       width=20, height=8, output.dir=output.dir)
      
      def.group.dotplot(enrich_gsea.result=gsea.all, showCategories=3,
                        method='GSEA: Reactome',
                        title=paste(selected.database, comparison.name,
                                    'GSEA_dotplot.pdf', sep = "_"),
                        width=6, height=4, output.dir=output.dir, func_colours=NULL,
                        groups=paste0(paste(val.1, collapse = '_'), ' vs ',
                                      paste(val.2, collapse = '_')))
      
      def.plot.barplot(enrich_gsea.result=gsea.all, showCategories=40,
                       method='GSEA: Reactome',
                       title=paste(selected.database, comparison.name,
                                   'GSEA_barplot.pdf', sep = "_"),
                       width=6, height=8, output.dir=output.dir, func_colours=NULL)
    }
  }
}

dataframe.pairs <- as.data.frame(df.pairs)

openxlsx::write.xlsx(dataframe.pairs, file.path(output.dir, 'Dendrogram_pairs.xlsx'))


