rm(list = ls(all=TRUE))
# Set seed before loading packages 
set.seed(1)

library('pathfindR')
library('clusterProfiler')
library('fgsea')
library('DOSE')
library(ggplot2)
library(xlsx)

library(ReactomePA)

def.run_pathfindR <- function(df_dge, custom.genes, custom.descriptions, pval.cut, 
                              multitest.method, output.dir) 
{
  custom_enrich_result <- pathfindR::run_pathfindR(
    df_dge, gene_sets = "Custom", custom_genes = custom.genes, 
    custom_descriptions = custom.descriptions, adj_method = multitest.method,
    p_val_threshold = pval.cut, # filter input dataframe
    enrichment_threshold = 0.05, # filter enrichment results 
    sig_gene_thr = 0.02, max_gset_size = Inf, # DO NOT LIMIT GENE SET SIZE
    output_dir = output.dir)
  
  ################### ---> convert gene ID to Symbol <--- ################### 
  custom_enrich_result = setreadable_pa(paenrich_object = custom_enrich_result)
  
  return(custom_enrich_result)
}


#' @param ranked.genes : all gene names of ranked geneList by e.g. log2FC
#' @param gene.sets : List of gene sets to check, e.g. Pathway names with associated genes
def.run_fgsea <- function(gene.sets, ranked.genes) 
{
  # GSEA with predefined pathways
  df_fgsea.res <- fgsea::fgsea(pathways = gene.sets, # List of gene sets to check.
                               stats    = ranked.genes, # named vector
                               minSize  = 10, maxSize  = 500, nPermSimple = 1000)
  
  return(df_fgsea.res)
}


#' @description performs a hypergeometric test comparing the set of "significant" genes against the "universe" (or background) genes
#' 
#' @param gene.names : all gene names of ranked geneList by e.g. log2FC
#' @param universe_genes : background genes / all measured genes in experiment
#' @param term.2.gene : Custom term (e.g. pathway names) to gene name (e.g. ENTREZID)
#' @param multitest.method : should be 'BH'
#' @param pval.cut : P-value cut-off parameter
def.run_enricher.ora <- function(gene.names, universe_genes, term.2.gene, multitest.method, pval.cut) 
{
  enricher.obj <- clusterProfiler::enricher(gene = gene.names,
                                            pvalueCutoff = pval.cut,
                                            pAdjustMethod = multitest.method,
                                            universe = universe_genes,
                                            minGSSize = 10,
                                            maxGSSize = 500,
                                            qvalueCutoff = 0.2,
                                            TERM2GENE = term.2.gene,
                                            TERM2NAME = NA)
  
  ################### ---> convert gene ID to Symbol <--- ################### 
  enricher.obj = setreadable_pa(paenrich_object = enricher.obj)
  
  return(enricher.obj)
}

def.run_enricher.gsea <- function(gene.list, term.2.gene, multitest.method, pval.cut) 
{
  enricher.obj <- clusterProfiler::GSEA(gene=gene.list,
                                        exponent = 1,
                                        minGSSize = 10,
                                        maxGSSize = 500,
                                        eps = 1e-10,
                                        pvalueCutoff = pval.cut,
                                        pAdjustMethod = multitest.method,
                                        TERM2GENE = term.2.gene,
                                        TERM2NAME = NA,
                                        verbose = TRUE,
                                        seed = TRUE,
                                        by = "fgsea",)
  
  ################### ---> convert gene ID to Symbol <--- ################### 
  enricher.obj = setreadable_pa(paenrich_object = enricher.obj)
  
  return(enricher.obj)
}


#' @param gene.names : all gene names of ranked geneList by e.g. log2FC
#' @param multitest.method : should be 'BH'
#' @param pval.cut : P-value cut-off parameter
#' @param ontology : one of "BP", "MF", and "CC" subontologies, or "ALL" for all three.
#' @param universe_genes : background genes / all measured genes in experiment
def.run.enrichgo.ora <- function(gene.names, universe_genes, multitest.method, pval.cut, ontology) 
{
  # Over-representation analysis
  enrichgo.obj <-clusterProfiler::enrichGO(
    gene = gene.names, OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID", ont=ontology,
    pvalueCutoff = pval.cut,  pAdjustMethod = multitest.method,
    universe = universe_genes, qvalueCutoff = 0.2, 
    minGSSize = 10, maxGSSize = 500,) 
  
  ################### ---> convert gene ID to Symbol <--- ################### 
  enrichgo.obj = setreadable_pa(paenrich_object = enrichgo.obj)
  
  return(enrichgo.obj)
}


#' @param gene.list : all genes ranked by e.g. log2FC
#' @param multitest.method : should be 'BH'
#' @param pval.cut : P-value cut-off parameter
#' @param ontology : one of "BP", "MF", and "CC" subontologies, or "ALL" for all three.
def.run.enrichgo.gsea <- function(gene.list, multitest.method, pval.cut, ontology) 
{
  # Gene set Enrichment Analysis
  enrichgo.obj <-clusterProfiler::gseGO(
    geneList = gene.list, OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID", ont=ontology,
    pvalueCutoff = pval.cut,  pAdjustMethod = multitest.method,
    minGSSize = 10, maxGSSize = 500,
    by = "fgsea") 
  
  ################### ---> convert gene ID to Symbol <--- ################### 
  enrichgo.obj = setreadable_pa(paenrich_object = enrichgo.obj)
  
  return(enrichgo.obj)
}


def.get_upregulated_genes <- function(func.degs) 
{
  func.degs <- func.degs[order(func.degs$sign.padj, decreasing = TRUE), ]
  de.diag <- func.degs$sign.padj
  names(de.diag) <- func.degs$entrezid
  de.diag.L_up <- de.diag[de.diag > -log10(0.05)]
  
  return(de.diag.L_up)
}

#' @author https://github.com/hms-dbmi/UpSetR/issues/85
fromList.v1 <- function (input) {
  # Same as original fromList()...
  elements <- unique(unlist(input))
  data <- unlist(lapply(input, function(x) {
    x <- as.vector(match(elements, x))
  }))
  data[is.na(data)] <- as.integer(0)
  data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(input), byrow = F))
  data <- data[which(rowSums(data) != 0), ]
  names(data) <- names(input)
  # ... Except now it conserves your original value names!
  row.names(data) <- elements
  return(data)
}

# ========================= Init Parameters ========================= #
output.dir <- '/Users/christina.hillig/R_studio/K_ERC_Grant/output'
output.dir <- file.path(output.dir, 'Pathways', 'Molecular_subtypes_res0.9_1vsRest', Sys.Date())
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

input.dir <- file.path('/Users/christina.hillig/R_studio/K_ERC_Grant/output', 'DGE_analysis', 'Molecular_subtypes_res0.9_1vsRest')


multitest.method <- 'BH'
pval.cut <- 0.05

date.data <- '2023-06-01' # '2023-02-24'

# ========================= Load Data ========================= #
df.dge.0 <- read.csv(file.path(input.dir, date.data, 'DEGs_E1_vs_Rest__Age_Sex_batchID_Molecular_subtype.csv'))
df.dge.1 <- read.csv(file.path(input.dir, date.data, 'DEGs_E2_vs_Rest__Age_Sex_batchID_Molecular_subtype.csv'))
df.dge.2 <- read.csv(file.path(input.dir, date.data, 'DEGs_E3_vs_Rest__Age_Sex_batchID_Molecular_subtype.csv'))

df.dge.3 <- read.csv(file.path(input.dir, date.data, 'DEGs_E4_vs_Rest__Age_Sex_batchID_Molecular_subtype.csv'))
df.dge.4 <- read.csv(file.path(input.dir, date.data, 'DEGs_E5_vs_Rest__Age_Sex_batchID_Molecular_subtype.csv'))
df.dge.5 <- read.csv(file.path(input.dir, date.data, 'DEGs_E6_vs_Rest__Age_Sex_batchID_Molecular_subtype.csv'))

df.dge.6 <- read.csv(file.path(input.dir, date.data, 'DEGs_E7_vs_Rest__Age_Sex_batchID_Molecular_subtype.csv'))
df.dge.7 <- read.csv(file.path(input.dir, date.data, 'DEGs_E8_vs_Rest__Age_Sex_batchID_Molecular_subtype.csv'))
df.dge.8 <- read.csv(file.path(input.dir, date.data, 'DEGs_E9_vs_Rest__Age_Sex_batchID_Molecular_subtype.csv'))

df.dge.9 <- read.csv(file.path(input.dir, date.data, 'DEGs_E10_vs_Rest__Age_Sex_batchID_Molecular_subtype.csv'))

df.dge.10 <- read.csv(file.path(input.dir, date.data, 'DEGs_E11_vs_Rest__Age_Sex_batchID_Molecular_subtype.csv'))
df.dge.11 <- read.csv(file.path(input.dir, date.data, 'DEGs_E12_vs_Rest__Age_Sex_batchID_Molecular_subtype.csv'))
df.dge.12 <- read.csv(file.path(input.dir, date.data, 'DEGs_E13_vs_Rest__Age_Sex_batchID_Molecular_subtype.csv'))



# ========================= Prepare Data ========================= #
# Add signed p.adj values
df.dge.0$sign.padj <- -sign(df.dge.0$log2FoldChange) * log10(df.dge.0$padj)
df.dge.1$sign.padj <- -sign(df.dge.1$log2FoldChange) * log10(df.dge.1$padj)
df.dge.2$sign.padj <- -sign(df.dge.2$log2FoldChange) * log10(df.dge.2$padj)

df.dge.3$sign.padj <- -sign(df.dge.3$log2FoldChange) * log10(df.dge.3$padj)
df.dge.4$sign.padj <- -sign(
  df.dge.4$log2FoldChange) * log10(df.dge.4$padj)
df.dge.5$sign.padj <- -sign(df.dge.5$log2FoldChange) * log10(df.dge.5$padj)

df.dge.6$sign.padj <- -sign(df.dge.6$log2FoldChange) * log10(df.dge.6$padj)
df.dge.7$sign.padj <- -sign(df.dge.7$log2FoldChange) * log10(df.dge.7$padj)
df.dge.8$sign.padj <- -sign(df.dge.8$log2FoldChange) * log10(df.dge.8$padj)

df.dge.9$sign.padj <- -sign(df.dge.9$log2FoldChange) * log10(df.dge.9$padj)
df.dge.10$sign.padj <- -sign(df.dge.10$log2FoldChange) * log10(df.dge.10$padj)
df.dge.11$sign.padj <- -sign(df.dge.11$log2FoldChange) * log10(df.dge.11$padj)
df.dge.12$sign.padj <- -sign(df.dge.12$log2FoldChange) * log10(df.dge.12$padj)

# Extract upregulated genes in lesion skin 
de.0.L_up <- def.get_upregulated_genes(func.degs = df.dge.0) 
de.1.L_up <- def.get_upregulated_genes(func.degs = df.dge.1) 
de.2.L_up <- def.get_upregulated_genes(func.degs = df.dge.2) 

de.3.L_up <- def.get_upregulated_genes(func.degs = df.dge.3) 
de.4.L_up <- def.get_upregulated_genes(func.degs = df.dge.4) 
de.5.L_up <- def.get_upregulated_genes(func.degs = df.dge.5) 


de.6.L_up <- def.get_upregulated_genes(func.degs = df.dge.6) 
de.7.L_up <- def.get_upregulated_genes(func.degs = df.dge.7) 
de.8.L_up <- def.get_upregulated_genes(func.degs = df.dge.8) 


de.9.L_up <- def.get_upregulated_genes(func.degs = df.dge.9) 
de.10.L_up <- def.get_upregulated_genes(func.degs = df.dge.10) 
de.11.L_up <- def.get_upregulated_genes(func.degs = df.dge.11) 
de.12.L_up <- def.get_upregulated_genes(func.degs = df.dge.12) 


# ========================= Perform Pathway analysis ========================= #
# ORA
pa.0.L_up <- ReactomePA::enrichPathway(
  gene=names(de.0.L_up), pvalueCutoff=pval.cut, readable=T, universe=df.dge.0$entrezid)
head(pa.0.L_up@result)
pa.1.L_up <- ReactomePA::enrichPathway(
  gene=names(de.1.L_up), pvalueCutoff=pval.cut, readable=T, universe=df.dge.1$entrezid)
pa.2.L_up <- ReactomePA::enrichPathway(
  gene=names(de.2.L_up), pvalueCutoff=pval.cut, readable=T, universe=df.dge.2$entrezid)


pa.3.L_up <- ReactomePA::enrichPathway(
  gene=names(de.3.L_up), pvalueCutoff=pval.cut, readable=T, universe=df.dge.3$entrezid)
pa.4.L_up <- ReactomePA::enrichPathway(
  gene=names(de.4.L_up), pvalueCutoff=pval.cut, readable=T, universe=df.dge.4$entrezid)
pa.5.L_up <- ReactomePA::enrichPathway(
  gene=names(de.5.L_up), pvalueCutoff=pval.cut, readable=T, universe=df.dge.5$entrezid)


pa.6.L_up <- ReactomePA::enrichPathway(
  gene=names(de.6.L_up), pvalueCutoff=pval.cut, readable=T, universe=df.dge.6$entrezid)
pa.7.L_up <- ReactomePA::enrichPathway(
  gene=names(de.7.L_up), pvalueCutoff=pval.cut, readable=T, universe=df.dge.7$entrezid)
pa.8.L_up <- ReactomePA::enrichPathway(
  gene=names(de.8.L_up), pvalueCutoff=pval.cut, readable=T, universe=df.dge.8$entrezid)

pa.9.L_up <- ReactomePA::enrichPathway(
  gene=names(de.9.L_up), pvalueCutoff=pval.cut, readable=T, universe=df.dge.9$entrezid)
pa.10.L_up <- ReactomePA::enrichPathway(
  gene=names(de.10.L_up), pvalueCutoff=pval.cut, readable=T, universe=df.dge.10$entrezid)
pa.11.L_up <- ReactomePA::enrichPathway(
  gene=names(de.11.L_up), pvalueCutoff=pval.cut, readable=T, universe=df.dge.11$entrezid)
pa.12.L_up <- ReactomePA::enrichPathway(
  gene=names(de.12.L_up), pvalueCutoff=pval.cut, readable=T, universe=df.dge.12$entrezid)


# Save results to .xlsx file
xlsx::write.xlsx(pa.0.L_up@result, file.path(output.dir, 'Cluster_E1.xlsx')) 
xlsx::write.xlsx(pa.1.L_up@result, file.path(output.dir, 'Cluster_E2.xlsx')) 
xlsx::write.xlsx(pa.2.L_up@result, file.path(output.dir, 'Cluster_E3.xlsx')) 
xlsx::write.xlsx(pa.3.L_up@result, file.path(output.dir, 'Cluster_E4.xlsx')) 
xlsx::write.xlsx(pa.4.L_up@result, file.path(output.dir, 'Cluster_E5.xlsx')) 
xlsx::write.xlsx(pa.5.L_up@result, file.path(output.dir, 'Cluster_E6.xlsx')) 
xlsx::write.xlsx(pa.6.L_up@result, file.path(output.dir, 'Cluster_E7.xlsx')) 
xlsx::write.xlsx(pa.7.L_up@result, file.path(output.dir, 'Cluster_E8.xlsx')) 
xlsx::write.xlsx(pa.8.L_up@result, file.path(output.dir, 'Cluster_E9.xlsx')) 
xlsx::write.xlsx(pa.9.L_up@result, file.path(output.dir, 'Cluster_E10.xlsx')) 
xlsx::write.xlsx(pa.10.L_up@result, file.path(output.dir, 'Cluster_E11.xlsx')) 
xlsx::write.xlsx(pa.11.L_up@result, file.path(output.dir, 'Cluster_E12.xlsx')) 
xlsx::write.xlsx(pa.12.L_up@result, file.path(output.dir, 'Cluster_E13.xlsx')) 


# ========================= Plot Pathway analysis ========================= #
library(clusterProfiler)
library(enrichplot)
library('org.Hs.eg.db')
# data(gcSample)
len <- 13
empty_list <- vector(mode = "list", length = len)
empty_list[[1]] <- names(de.0.L_up)
empty_list[[2]] <- names(de.1.L_up)
empty_list[[3]] <- names(de.2.L_up)

empty_list[[4]] <- names(de.3.L_up)
empty_list[[5]] <- names(de.4.L_up)
empty_list[[6]] <- names(de.5.L_up)

empty_list[[7]] <- names(de.6.L_up)
empty_list[[8]] <- names(de.7.L_up)
empty_list[[9]] <- names(de.8.L_up)

empty_list[[10]] <- names(de.9.L_up)
empty_list[[11]] <- names(de.10.L_up)
empty_list[[12]] <- names(de.11.L_up)
empty_list[[13]] <- names(de.12.L_up)

names(empty_list) <- c('E1', 'E2', 'E3', 'E4', 'E5', 'E6', 'E7', 'E8', 'E9', 'E10', 'E11', 'E12', 'E13')
# names(empty_list) <- as.integer(names(empty_list)) + 1
# xx <- clusterProfiler::compareCluster(
#   empty_list, fun="enrichKEGG", organism="human", pvalueCutoff=0.05, universe=df.dge.12$entrezid)
# xx <- enrichplot::pairwise_termsim(xx) 

xx.pa <- clusterProfiler::compareCluster(
  empty_list, fun="enrichPathway", pvalueCutoff=0.05, universe=df.dge.12$entrezid)
xx.pa <- enrichplot::pairwise_termsim(xx.pa)   

xx.go <- clusterProfiler::compareCluster(
  empty_list, fun="enrichGO", pvalueCutoff=0.05, OrgDb='org.Hs.eg.db', universe=df.dge.12$entrezid)
xx.go <- enrichplot::pairwise_termsim(xx.go)   


# Save results to .xlsx file
xlsx::write.xlsx(xx.pa@compareClusterResult, file.path(output.dir, 'Clusterd_Pathways.xlsx')) 
xlsx::write.xlsx(xx.go@compareClusterResult, file.path(output.dir, 'Clusterd_Goterms.xlsx')) 
# xlsx::write.xlsx(xx@compareClusterResult, file.path(output.dir, 'Clusterd_Pathways_KEGG.xlsx')) 


# colors = c(palette.colors(palette = "Okabe-Ito"), 
#            palette.colors(palette = "Set2"))[1:length(xx@geneClusters)]
colors <- c("#006ddb", "#b6dbff", "#004949", "#009292", 
            "#ff6db6", "#490092", "#b66dff", "#000000", "#920000", 
            "#E69F00", "#D55E00", "#8B4513", "#999999", 
            "#6db6ff", "#924900", "#db6d00", "#ffff6d", "#24ff24", 
            "#ffb6db")[1:length(xx.pa@geneClusters)]
# colors <- c("#000000","#004949","#009292","#ff6db6", "#ffb6db",
#             "#490092","#006ddb","#b66dff","#6db6ff","#b6dbff",
#             "#920000","#924900","#db6d00", "#24ff24","#ffff6d")[1:length(xx.pa@geneClusters)]
# p1 <- enrichplot::emapplot(
#   xx, color = 'p.adjust',
#                showCategory = 5,
#                node_scale = 0.1,
#                line_scale = 0.5,
#                min_edge = 0.1,
#                cex_label_category = 0.5,
#                cex_category = 1,
#                cex_line = 0.5) + 
#   scale_fill_manual(values=unname(colors))
# p1$labels$fill <- "Endotypes"
# 
# pdf(file = file.path(output.dir, paste0('Clustered_KEGG.pdf')), 
#     width = 6, height = 6)
# print(p1)
# dev.off()

# Pathway EMAP plot
p2 <- emapplot(xx.pa,
               showCategory = 4,
               node_scale = 0.1,
               line_scale = 0.01,
               min_edge = 0.1,
               node_label_size = 0.1,
               cex_label_category = 1.3,  # Label size
               repel=FALSE, 
               direction='both', # "both", "x", or "y" – direction in which to adjust position of labels
               cex_category = 1,
               cex_line = 0.8, # Scale of line width
               label_style='shadowtext', 
               shadowtext=FALSE, 
               nWords = 3, 
               cex_pie2axis = 1, 
               legend_n = 4)  +  
  scale_fill_manual(values=unname(colors))  + 
  ggplot2::theme(legend.text = element_text(size = 10), legend.title = element_text(size = 12))
p2$labels$fill <- "Endotypes"

pdf(file = file.path(output.dir, paste0('Clustered_Pathway.pdf')), 
    width = 10, height = 7)
print(p2)
dev.off()


p3 <- emapplot(xx.go,
               showCategory = 5,
               node_scale = 0.1,
               line_scale = 0.5,
               min_edge = 0.1,
               cex_label_category = 0.5,
               cex_category = 1,
               cex_line = 0.5) +  ggplot2::theme(legend.text = element_text(size = 8), 
                                                 legend.title = element_text(size = 10)) + 
  scale_fill_manual(values=unname(colors))
p3$labels$fill <- "Endotypes"
pdf(file = file.path(output.dir, paste0('Clustered_Goterms.pdf')), 
    width = 6, height = 6)
print(p3)
dev.off()

# Plot pathways which are unique to each cluster
len <- 13
pas_list <- vector(mode = "list", length = len)
pas_list[[1]] <- pa.0.L_up@result$Description
pas_list[[2]] <- pa.1.L_up@result$Description
pas_list[[3]] <- pa.2.L_up@result$Description

pas_list[[4]] <- pa.3.L_up@result$Description
pas_list[[5]] <- pa.4.L_up@result$Description
pas_list[[6]] <- pa.5.L_up@result$Description

pas_list[[7]] <- pa.6.L_up@result$Description
pas_list[[8]] <- pa.7.L_up@result$Description
pas_list[[9]] <- pa.8.L_up@result$Description

pas_list[[10]] <- pa.9.L_up@result$Description
pas_list[[11]] <- pa.10.L_up@result$Description
pas_list[[12]] <- pa.11.L_up@result$Description
pas_list[[13]] <- pa.12.L_up@result$Description
names(pas_list) <- names(empty_list)
df.matrix <- fromList.v1(pas_list)

# In order to find unique genes rowSums have to be equal to 1
# v.tmp.rowsums <- venn.up_down[rowSums(venn.up_down) ==1, ]
# However, its tricky to find unique pathways therefore we include pathways 
# which are shared between max. 6 comparisons
v.tmp.rowsums <- df.matrix[(rowSums(df.matrix) <= 6) & (rowSums(df.matrix) > 0), ]
signautre.pas = hash::hash()
signautre.pas.vector = c()
for (cluster.names in colnames(df.matrix))
{
  # Select genes from condition of interest, 1 means gene is DEx in L vs NL comparison
  cluster_sig_pas <- row.names(
    v.tmp.rowsums[v.tmp.rowsums[[cluster.names]] == 1, ])
  
  signautre.pas[[cluster.names]] <- cluster_sig_pas
  signautre.pas.vector <- c(signautre.pas.vector, cluster_sig_pas)
  
}
signautre.pas.vector <- unique(signautre.pas.vector)



#  Does nit work atm
# pdf(file = file.path(output.dir, paste0('Clustered_Pathway_treeplot.pdf')), 
#     width = 16, height = 10)
# print(
#   treeplot(xx.pa, showCategory = 8, hclust_method = "ward.D", nCluster = 5, nWords = 3,
#            label_format = 0.5, cex_category = 1, fontsize=4, offset_tiplab = 0.7,
#            offset = 1, extend = 0.3, hexpand =0.05)
# )
# dev.off()
# 
# 
# pdf(file = file.path(output.dir, paste0('Clustered_Goterms_treeplot.pdf')), 
#     width = 14, height = 10)
# print(
#   treeplot(xx.go, showCategory = 8, hclust_method = "ward.D", nCluster = 5, nWords = 3,
#            label_format = 0.5, cex_category = 1, fontsize=4, offset_tiplab = 0.7,
#            offset =1, extend = 0.3, hexpand =0.05, align='both') +
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# )
# dev.off()





