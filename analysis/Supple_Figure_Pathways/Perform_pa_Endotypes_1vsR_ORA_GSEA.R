# Script to perform DGE analysis between branches of a node
rm(list = ls(all=TRUE))

set.seed(1)

library('xlsx')
library(openxlsx)
library(data.table)
library("dendextend")
# Pathway analysis libraries
library(ReactomePA)
library('fgsea')
library('DOSE')
library('org.Hs.eg.db')

general_path = '/Volumes/CH__data/Projects/Eyerich_AG_projects/BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer/analysis/Molecular_subtypes'
source(file.path(general_path, 'r_scripts/Pathway_analysis/ORA_GSEA_plots.R'))

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
selected.database <- 'ReactomePA'
fdr.value <- 0.05
padj.cut_degs <- 0.05
l2fc.factor <- 1
minGSSize <- 10

multitest.method <- 'BH'
pval.cut <- 0.05

date.data <- '2023-06-01'
general.dir <- file.path(general_path, 'output')
input.dir <- file.path(general.dir, 'DGE_analysis', 'Molecular_subtypes_res0.9_1vsRest', date.data)
output.dir <- file.path(general.dir, 'Pathways', 'Endotypes_1_vs_R_ORA_GSEA', Sys.Date())
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)


# Use list.files to find all files in subdirectories
all_files <- list.files(input.dir, recursive = TRUE, full.names = TRUE)

# Filter the list to include only CSV files
csv_files <- all_files[grep("\\.csv$", all_files)]

hash.dict <- hash::hash()
# Save upregulated gene lists to list
DEGs_entrezID <- vector(mode = "list", length = length(csv_files))
DEGs_hgnc <- vector(mode = "list", length = length(csv_files))
DEGs.ranked_genes <- vector(mode = "list", length = length(csv_files))
ora_list <- vector(mode = "list", length = length(csv_files))
gsea_list <- vector(mode = "list", length = length(csv_files))
obj.gsea_list <- vector(mode = "list", length = length(csv_files))

message('Pathway enrichment analysis of L skin biopsies of one endotype vs. all others')
for (file_name in csv_files) 
{
  df.dge <- read.csv(file = file_name, row.names = 1)
  comparison.name <- strsplit(strsplit(file_name, split='/')[[1]][13], split='__')[[1]][1]
  val.1 <- strsplit(comparison.name, split = '_')[[1]][2]
  val.2 <- strsplit(comparison.name, split = '_')[[1]][4]
  
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
    gene=names(ranked_genes.gs_enID.up[[2]]), pvalueCutoff=pval.cut, 
    readable=T, universe = df.dge$entrezid)
  ora_list[[paste(val.1, collapse = '_')]] <- pa.L_up@result[
    pa.L_up@result$p.adjust < fdr.value, ]$Description
  
  pa.L_down <- ReactomePA::enrichPathway(
    gene=names(ranked_genes.gs_enID.down[[2]]), pvalueCutoff=pval.cut, 
    readable=T, universe = df.dge$entrezid)
  ora_list[[paste(val.2, collapse = '_')]] <- pa.L_down@result[
    pa.L_down@result$p.adjust < fdr.value, ]$Description
  
  # GSEA
  gsea.all <- ReactomePA::gsePathway(
    ranked_genes.enID.all, organism = "human", eps=0,
    pvalueCutoff=pval.cut, pAdjustMethod=multitest.method, verbose=TRUE,
    minGSSize = minGSSize, by = "fgsea", seed=TRUE)
  
  gsea.all <- DOSE::setReadable(gsea.all, 'org.Hs.eg.db', 'ENTREZID')
  gsea_list[[comparison.name]] <- gsea.all@result[
    gsea.all@result$p.adjust < pval.cut, ]$Description
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
# Extract all pathways which are upregulated in endotypes
ora_list[[paste(val.1, collapse = '_')]]
# Positive enrichement score 
gsea_list[[comparison.name]] 
gsea_list.v2 <- gsea_list[names(gsea_list)[names(gsea_list) != '']]
gsea.table.up <- fromList.v1(gsea_list.v2)

signature_gsea_pathways <- vector(mode = "list")
for (comparison.name in names(gsea_list)[names(gsea_list) != '']) 
{
  gsea.tmp.rowsums <- gsea.table.up[rowSums(gsea.table.up) == 1, ]
  condition_sig_updown <- row.names(
    gsea.tmp.rowsums[gsea.tmp.rowsums[[comparison.name]] == 1, ])
  
  signature_gsea_pathways[[comparison.name]] <- condition_sig_updown
  
  xlsx::write.xlsx(condition_sig_updown, file = file.path(
    output.dir,  paste('Unique_GSEA', '.xlsx', sep='')), 
    sheetName=comparison.name, append=TRUE)
}

# TODO highlighy the composition
plot1 <- UpSetR::upset(gsea.table.up, nsets=13, keep.order = F,
                       nintersects=10, order.by = "freq", decreasing = TRUE)
plot2 <- UpSetR::upset(gsea.table.up, nsets=13, keep.order = T, 
                       nintersects=12, set.colors=colors)

# For ORA results:
ora.up.filtered <- ora_list[names(ora_list)[names(ora_list) != "" & names(ora_list) != 'Rest']]
ora.table.up <- fromList.v1(ora.up.filtered)

pdf(file = file.path(output.dir, paste0('Upsetplot_unique_ORA_pathways.pdf')), 
    width = 6, height = 6)
print(UpSetR::upset(ora.table.up, nsets=13, keep.order = T, nintersects=12))
dev.off()

signature_ora_pathways <- vector(mode = "list")
for (comparison.name in names(ora.up.filtered)) 
{
  ora.tmp.rowsums <- ora.table.up[rowSums(ora.table.up) == 1, ]
  condition_sig_updown <- row.names(
    ora.tmp.rowsums[ora.tmp.rowsums[[comparison.name]] == 1, ])
  
  signature_ora_pathways[[comparison.name]] <- condition_sig_updown
  
  xlsx::write.xlsx(condition_sig_updown, file = file.path(
    output.dir,  paste('Unique_ORA', '.xlsx', sep='')), 
    sheetName=comparison.name, append=TRUE)
}

# reorder
reordered.DEGs_entrezID <- DEGs_entrezID[names(DEGs_entrezID)[names(DEGs_entrezID) != "" & names(
  DEGs_entrezID) != 'Rest']][c('E1', 'E2', 'E3', 'E4', 'E5', 'E6', 'E7', 'E8', 
                               'E9', 'E10', 'E11', 'E12', 'E13')]

xx.pa <- clusterProfiler::compareCluster(
  reordered.DEGs_entrezID, fun="enrichPathway", pvalueCutoff=0.05, 
  universe=df.dge$entrezid, readable=T)
# xx.pa <- clusterProfiler::compareCluster(
#   DEGs.ranked_genes, fun="gsePathway", pvalueCutoff=pval.cut, 
#   pAdjustMethod=multitest.method, organism = "human", eps=0, 
#   minGSSize = minGSSize, by = "fgsea", seed=TRUE)

xx.pa@compareClusterResult$Cluster <- factor(
  xx.pa@compareClusterResult$Cluster, 
  levels=c('E1', 'E2', 'E3', 'E4', 'E5', 'E6', 'E7', 'E8', 'E9', 'E10', 'E11', 'E12', 'E13'))
xx.pa <- enrichplot::pairwise_termsim(xx.pa)   



colors <- c("#006ddb", "#b6dbff", "#004949", "#009292", 
            "#ff6db6", "#490092", "#b66dff", "#000000", "#920000",  
            "#E69F00", "#D55E00", "#8B4513", "#999999")

names(colors) <- c('E1', 'E2', 'E3', 'E4', 'E5', 'E6', 'E7', 'E8', 
                   'E9', 'E10', 'E11', 'E12', 'E13')

xx.pa@compareClusterResult$Cluster <- factor(
  xx.pa@compareClusterResult$Cluster, 
  levels=c('E1', 'E2', 'E3', 'E4', 'E5', 'E6', 'E7', 'E8', 'E9', 'E10', 'E11', 'E12', 'E13'))
p2 <- emapplot(xx.pa,
               showCategory = 5,
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
               legend_n = 4)   + 
  ggplot2::theme(legend.text = element_text(size = 10), legend.title = element_text(size = 12)) +
  # ggplot2::scale_fill_discrete(limits = c('E1', 'E2', 'E3', 'E4', 'E5', 'E6', 'E7', 'E8', 
  #                                         'E9', 'E10', 'E11', 'E12', 'E13')) +
  scale_fill_manual(values=colors) 
p2$labels$fill <- "Endotypes"

pdf(file = file.path(output.dir, paste0('Clustered_Pathway.pdf')), 
    width = 6, height = 6)
print(p2)
dev.off()

# ------------- Pathway Family
pathways.id_name <- toTable(reactome.db::reactomePATHNAME2ID)
pathways.id_name <- pathways.id_name[
  grep("Homo sapiens: ", iconv(pathways.id_name$path_name)), ]
pathways.id_name <- pathways.id_name[grep("-HSA-", iconv(pathways.id_name$DB_ID)), ]
pathways.id_name <- pathways.id_name %>% dplyr::rename(description = path_name)

pathways.id_entrezID <- toTable(reactome.db::reactomePATHID2EXTID)
pathways.id_entrezID <- pathways.id_entrezID[grep("-HSA-", iconv(pathways.id_entrezID$DB_ID)), ]
pathways.id_entrezID <- pathways.id_entrezID %>% dplyr::rename(entrezid = gene_id)

# annotated with the SYMBOL and ENSEMBL identifiers associated with each Entrez id
hsa_reactome_anno <- pathways.id_entrezID %>%
  dplyr::mutate(
    symbol = mapIds(org.Hs.eg.db, entrezid, "SYMBOL", "ENTREZID"),
    ensembl = mapIds(org.Hs.eg.db, entrezid, "ENSEMBL", "ENTREZID")
  )

# join dataframes
ractome_pas <- dplyr::left_join(hsa_reactome_anno, pathways.id_name)

ractome_pas$description <- gsub("Homo sapiens: ", "", as.character(ractome_pas$description))

# get hierarchy of pathways
if (file.exists('/Users/christina.hillig/Downloads/ReactomePathwaysRelation_Description.csv')) 
{
  df.hierarchy <- read.csv('/Users/christina.hillig/Downloads/ReactomePathwaysRelation_Description.csv')
} else {
  df.hierarchy <- read.table('/Users/christina.hillig/Downloads/ReactomePathwaysRelation.txt')
  
  high.level.description <- c()
  for (pa.events in unique(df.hierarchy$V1)) 
  {
    tmp <- ReactomeContentService4R::query(id = pa.events)
    high.level.description <- c(high.level.description, tmp$displayName)
  }
  
  df.hierarchy$description <- 'Unknown'
  counter <- 1
  for (pa.events in unique(df.hierarchy$V1)) {
    mask <- df.hierarchy$V1 == pa.events
    df.hierarchy$description[mask] <- rep(high.level.description[counter], sum(mask))
    counter <- counter + 1
  }
  write.csv(df.hierarchy, file='/Users/christina.hillig/Downloads/ReactomePathwaysRelation_Description.csv')
}


# Get frequency on how many high level pathways have been found by low level pathways
# Only significant pathways
sig.ids <- xx.pa@compareClusterResult[xx.pa@compareClusterResult$p.adjust < 0.0001, ]
freq.pathways <- c()
mask.pas <- c()
for (counter in 1:dim(sig.ids)[1]) {
  pa.tmp <- df.hierarchy$description[df.hierarchy$V2 == sig.ids$ID[counter]]
  freq.pathways <- c(freq.pathways, pa.tmp)
  mask.pas <- c(mask.pas, sig.ids$ID[counter] %in% df.hierarchy$V2)
}

df.pa.clusters <- data.frame(Cluster=sig.ids$Cluster[mask.pas], ID=freq.pathways)
# Create a dataframe of the form from (cluster), to (pathway), values (frequency)
# Group by Cluster and ID, then count the occurrences
result <- df.pa.clusters %>%
  dplyr::group_by(Cluster, ID) %>%
  dplyr::count() %>%
  dplyr::ungroup() %>%
  dplyr::arrange(ID)  # Optional: Sort by ID

# Rename columns
names(result) <- c("from", "to", "value")
result$from <- factor(result$from, levels=c('E1', 'E2', 'E3', 'E4', 'E5', 'E6', 'E7', 'E8', 'E9', 'E10', 'E11', 'E12', 'E13'))
grid.col = c(
  E1 = "#006ddb", E2 = "#b6dbff", E3 = "#004949", E4 = "#009292", E5 = "#ff6db6", E6 = "#490092", 
  E7="#b66dff", E8="#000000", E9="#920000", E10="#E69F00", E11="#D55E00", E12="#8B4513", E13="#999999" )

library(circlize)

pdf(file = file.path(output.dir, paste0('ChordDiaggram_Pathways.pdf')), 
    width = 8, height = 8)
chordDiagramFromDataFrame(
  result, grid.col=grid.col, annotationTrack = "grid",
  order=union(c('E1', 'E2', 'E3', 'E4', 'E5', 'E6', 'E7', 'E8', 
                'E9', 'E10', 'E11', 'E12', 'E13'), 
              result$to))
  

# first = tapply(result$from, result$to, function(x) x[1])
# last = tapply(result$from, result$to, function(x) x[length(x)])
# for(i in seq_along(first)) {
#   start.degree = get.cell.meta.data("cell.start.degree", sector.index = first[i], track.index = 1)
#   end.degree = get.cell.meta.data("cell.end.degree", sector.index = last[i], track.index = 1)
#   rou1 = get.cell.meta.data("cell.bottom.radius", sector.index = first[i], track.index = 1)
#   rou2 = get.cell.meta.data("cell.top.radius", sector.index = last[i], track.index = 1)
#   draw.sector(start.degree, end.degree, rou1, rou2, border = NA, col = "white")
# }

for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
  circos.text(mean(xlim), ylim[1], si, facing = "clockwise", adj = c(0, 0.5),
              niceFacing = TRUE, cex = 1.25, col = "black", 
              sector.index = si, track.index = 1)
}

# highlight.sector(list('E11', 'E12', 'E13'), track.index = 1, col = "#6db6ff", 
#                  text = "A", cex = 0.8, text.col = "black", niceFacing = TRUE)
# highlight.sector(list('E1', 'E2', 'E3', 'E4'), track.index = 1, col = "darkorange", 
#                  text = "A", cex = 0.8, text.col = "black", niceFacing = TRUE)
# highlight.sector(list('E5', 'E6', 'E7', 'E8', 'E9', 'E10'), track.index = 1, col = "#920000", 
#                  text = "A", cex = 0.8, text.col = "black", niceFacing = TRUE)

dev.off()


                 




# 
# 
# library(GOplot)
# 
# df.toy <- df.dge
# df.toy$ID <- df.toy$hgnc_symbol
# df.toy$log2FC <- df.toy$log2FoldChange
# df.toy <- df.toy[(df.toy$log2FoldChange > 1) & (df.toy$padj < 0.05), ]
# 
# enrichgo.obj <- clusterProfiler::gseGO(
#   gene = ranked_genes.enID.all,  OrgDb = org.Hs.eg.db,
#   keyType = "ENTREZID", ont='ALL',
#   pvalueCutoff = pval.cut,  pAdjustMethod = multitest.method,
#   minGSSize = 10, maxGSSize = 500,
#   by = "fgsea") 
# 
# enrichgo.obj@result$Term <- pa.L_up@result$Description
# enrichgo.obj@result$adj_pval <- pa.L_up@result$p.adjust  
# enrichgo.obj@result$Genes <- pa.L_up@result$core_enrichment  
# 
# circ <- circle_dat(enrichgo.obj, df.dge)
# 
# # Now it is time to generate the binary matrix
# chord <- chord_dat(circ, df.toy, pa.L_up@result$Description[pa.L_up@result$p.adjust < 0.05])



