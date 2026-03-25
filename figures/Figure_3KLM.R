# Metabolism analysis
# Create barplots and Cnetplots for each file

rm(list = ls())
set.seed(1)

library('xlsx')
library('dplyr')
library('hash')
library('multienrichjam')

# Plotting libraries
library('ggplot2')
library('ggnewscale')
library('enrichplot')


general_path = '/Volumes/CH__data/Projects/Eyerich_AG_projects/BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer/analysis'
# 1. Create output directory
save_dir <- file.path(general_path, 'output', 'Figure_3KLM', Sys.Date())
dir.create(save_dir, showWarnings = FALSE,  recursive = TRUE)

# 2. Load .csv files
input.dir.pas <- file.path(general_path, 'output/Pathways', 
                           'Endotypes_Dendrogram/2023-06-22')
input.dir.degs <- file.path(general_path, 'output/DGE_analysis', 
                            'Endotypes_Dendrogram/2023-08-29')

# Get all files in subfolders
all_files.pas = c("GSEA__E5_E6_E7_E8_E9_E10__ vs __E11_E12_E13_E1_E2_E3_E4__Age_Sex_batchID_subtypes.xlsx",
                  "GSEA__E11_E12_E13__ vs __E1_E2_E3_E4__Age_Sex_batchID_subtypes.xlsx",
                  "GSEA__E1_E2__ vs __E3_E4__Age_Sex_batchID_subtypes.xlsx")
all_files.degs = c("DEGs_E5_E6_E7_E8_E9_E10__ vs __E11_E12_E13_E1_E2_E3_E4__Age_Sex_batchID_subtypes.xlsx",
                   "DEGs_E11_E12_E13__ vs __E1_E2_E3_E4__Age_Sex_batchID_subtypes.xlsx",
                   "DEGs_E1_E2__ vs __E3_E4__Age_Sex_batchID_subtypes.xlsx" )

# read out all corresponding deg file names of pas
for (file in all_files.pas) 
{
  group_one <- strsplit(file, '__')[[1]][2]
  group_two <- strsplit(file, '__')[[1]][4]
  # why does this work!?
  temp.deg_filename <- all_files.degs[grep(paste0(group_one, '.*', group_two), all_files.degs)]
  
  # create per file barplot or dotplot
  df.degs_US <- xlsx::read.xlsx(file.path(
    input.dir.degs, temp.deg_filename), 
    sheetIndex = 1, header = TRUE, 
    Comparision=paste0(group_one, '_vs_', group_two), )
  
  df.all__temp <- xlsx::read.xlsx(file.path(
    input.dir.pas, file), 
    sheetIndex = 1, header = TRUE, 
    Comparision=paste0(group_one, '_vs_', group_two), )
  df.all__temp$NA. <- NULL
  
  # df.all__temp <- def.plot_pas(file.dir=input.dir.pas, func_stimulus=stimuli)
  
  # Read out all genes in these pathways and store in a list
  gene.lists <- strsplit(df.all__temp$core_enrichment, "/")
  gene.list <- unique(unlist(gene.lists, recursive = FALSE))
  
  # Get log2FC, p-value and p.adj value for gene in pahtways of interest
  ind.genes.deg <- match(gene.list, df.degs_US$hgnc_symbol)
  df.degs_goi <- df.degs_US[ind.genes.deg, ]
  
  # Save as lists separately
  write.xlsx(df.degs_goi, file.path(save_dir, paste0('Genes__GSEA__', temp.deg_filename)))
  
  # Convert GeneRatio to numeric
  df.all__temp$NES <- as.numeric(sapply(df.all__temp$NES, function(x) eval(parse(text=x))))
  
  # Create cnet- and barplot
  # Barplot
  png(
    filename = file.path(
      save_dir,
      paste('Barplot', paste0(strsplit(temp.deg_filename, '[.]')[[1]][1], "__GSEA.png"))),
    width = 8, height = 8, units = "in", bg = "white", res = 300)
  p1 <- ggplot(df.all__temp, aes(x = Description, y = NES, fill = p.adjust)) +
    geom_bar(stat = "identity", position = "dodge", width=0.5) + coord_flip() +
    theme(text = element_text(size=14)) # , aspect.ratio = 2/1
  print(p1)
  dev.off()
  
  # Cnetplot
  if (length(df.all__temp$Description) > 2) 
  {
    # Convert df to enrichresult
    enrich.temp <- multienrichjam::enrichDF2enrichResult(df.all__temp,
                                                         pvalueCutoff = 0.05,
                                                         pAdjustMethod = "BH",
                                                         keyColname = "ID",
                                                         geneColname = 'core_enrichment',
                                                         geneHits = "setSize",
                                                         geneRatioColname = 'NES',
                                                         geneDelim = "[,/ ]+",
                                                         pvalueColname = 'p.adjust',
                                                         descriptionColname = "Description",
                                                         msigdbGmtT = NULL,
                                                         verbose = FALSE,)
    # create named vector
    ranked.genes <- df.degs_goi$log2FoldChange
    names(ranked.genes) <-   df.degs_goi$hgnc_symbol
    ranked.genes <- ranked.genes[order(factor(names(ranked.genes), levels=enrich.temp@gene))]
    
    png(
      filename = file.path(
        save_dir,
        paste('Cnetplot', paste0(strsplit(temp.deg_filename, '[.]')[[1]][1], "__GSEA.png"))),
      width = 16, height = 8, units = "in", bg = "white", res = 300)
    p2 <- enrichplot::cnetplot(enrich.temp, 
                               showCategory = length(enrich.temp@result$Description), 
                               foldChange = ranked.genes, layout = "kk")
    print(p2)
    dev.off() 
  }
  
}


