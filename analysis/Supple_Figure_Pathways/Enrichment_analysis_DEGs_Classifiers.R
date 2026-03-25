# Figure 8
rm(list = ls())
set.seed(42)  # important: as permutations and random initialization are used
options(java.parameters = "-Xmx4g")
library(fgsea)
# library(ggplot2)
library('xlsx')
library(stringr)

# 1. Load patient Bulk DEGs of:
# - Pso, L vs Pso NL
# - LP & LE, L vs LP & LE, NL
# - AE, L vs AE, NL
# 2. Rank List by log2FC or signed p.adj value
# 3. Gene Set:
# - load DEGs from 2D skin models
# - filter DEGs by up/down regulated genes
# - create df with column names up / down for each disease
# 4. Perform fgsea

path.general <- '/Volumes/CH__data/Projects/Eyerich_AG_projects/BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer/analysis/Molecular_subtypes'


for (classifier in c('IL23_IL17', 'TNF'))
{
  # 1. Create output directory
  save_dir <- file.path(
    path.general, 'output', "Enrichment_analysis", classifier, Sys.Date())
  dir.create(save_dir, showWarnings = FALSE,  recursive = TRUE)
  
  if (classifier == 'TNF') 
  {
    df <- read.xlsx('/Volumes/CH__data/Projects/Eyerich_AG_projects/BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer/analysis/Molecular_subtypes/output/Compare_GEx_TNF_models/2024-04-11/Feature_importance.xlsx',
                    sheetIndex = 1)
    df.degs <- readxl::read_excel(
      '/Volumes/CH__data/Projects/Eyerich_AG_projects/BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer/analysis/Molecular_subtypes/output/DGE_analysis/Endotypes_Dendrogram/2023-08-29/DEGs_E12__ vs __E13__Age_Sex_batchID_subtypes.xlsx', 
      sheet = 1)
    df.degs$gene_symbol <- df.degs$hgnc_symbol
    
  } else 
  {
    # df <- read.csv('/Volumes/CH__data/Projects/Eyerich_AG_projects/BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer/analysis/Molecular_subtypes/output/Compare_GEx_IL23_IL17_models/2023-12-28/UpSetPlot_features_modelnames.csv')
    df <- read.xlsx('/Volumes/CH__data/Projects/Eyerich_AG_projects/BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer/analysis/Molecular_subtypes/output/Compare_GEx_IL23_IL17_models/2024-04-11/Feature_importance.xlsx', 
                    sheetIndex = 1)
    df.degs <- readxl::read_excel('/Volumes/CH__data/Thesis/writing/figures/Chapter_4/files/DGE_results.xlsx',
                                  sheet = 'E12_E13vs.E8_E11')
    df.degs$gene_symbol <- df.degs$Gene_name
    # Drop duplicated genes
    df.degs <- df.degs[!duplicated(df.degs$gene_symbol), ]
  }
  
  # features <- gsub("\\[|\\]", "", df$features)
  # features <- gsub("'", "", features)
  # named_list.features <- stringr::str_split(features, ', ')
  # 
  # names(named_list.features) <- df$model
  
  message(paste0('DEGs: ', classifier))
  # add singed p-value column
  df.degs$signed.padj <- sign(df.degs$log2FoldChange) * df.degs$padj
  # 1. Rank List
  df.ranked_degs <- df.degs[order(df.degs$signed.padj), ]
  # names vector
  ranked.genes <- df.ranked_degs$log2FoldChange
  names(ranked.genes) <- df.ranked_degs$gene_symbol
  
  
  message('Start GSEA analysis of DEGs using as gene sets of the classifier \n')
  for (val in unique(df$Model))
  {
    message('===========================')
    # 2. Gene Set
    message(paste0('Gene set: ', val))
    # list.classifier_genes <- named_list.features[[val]]
    # genesets <- list(list.classifier_genes)
    # names(genesets) <- c(val)
    
    df_tmp <- df[df$Model == val, ]
    
    # 3.1 Get up/down regulated genes
    class_1_features <- df_tmp[df_tmp$direction > 0, ]$features
    class_0_features <- df_tmp[df_tmp$direction  < 0, ]$features
    # 3.2 Create new named list
    genesets <- list(class_1_features, class_0_features)
    names(genesets) <- c(paste0(val, "_class_1"),
                         paste0(val, "_class_0"))
    
    # 4. fgsea
    # Reactome: pathways <- reactomePathways(names(exampleRanks))
    # pathways <- gmtPathways(gmt.file)
    fgseaRes <- fgsea::fgsea(pathways = genesets, # Classifier gene sets
                             stats    = ranked.genes, # DEGs
                             minSize  = 1, maxSize  = Inf, nPermSimple = 1000)
    
    # A table with GSEA results.
    # Each row corresponds to a tested pathway. The columns are the following
    # pathway – name of the pathway as in 'names(pathway)';
    # pval – an enrichment p-value;
    # padj – a BH-adjusted p-value;
    # log2err – the expected error for the standard deviation of the P-value logarithm.
    # ES – enrichment score, same as in Broad GSEA implementation;
    # NES – enrichment score normalized to mean enrichment of random samples of the same size;
    # size – size of the pathway after removing genes not present in 'names(stats)'.
    # leadingEdge – vector with indexes of leading edge genes that drive the enrichment, see http://software.broadinstitute.org/gsea/doc/GSEAUserGuideTEXT.htm#_Running_a_Leading.
    
    pdf(file = file.path(
      save_dir,
      paste('Enrichmentplot', paste0(val, '__class_1', ".pdf"))),
      width = 4, height =4)
    p.up <- plotEnrichment(genesets[[paste0(val, "_class_1")]], ranked.genes) +
      ggplot2::labs(title= paste0("Predictive features of class 1 in DEGs")) +
      ggplot2::annotate(
        "label", x=2000, y=-0.6, label.size = NA,
        label=paste0("NES: ", sprintf(fgseaRes$NES[1], fmt = '%#.2f'), "\npadj: ", 
                     sprintf(fgseaRes$padj[1], fmt = '%#.2e')))
    
    print(p.up)
    dev.off()
    
    pdf(file = file.path(
      save_dir,
      paste('Enrichmentplot', paste0(val, '__class_0', ".pdf"))),
      width = 4, height =4)
    p.down <- 
      plotEnrichment(genesets[[paste0(val, "_class_0")]], ranked.genes) +
      ggplot2::labs(title=paste0("Predictive features of class 0 in DEGs")) +
      ggplot2::annotate(
        "label", x=12000, y=0.6, label.size = NA,
        label=paste0("NES: ", sprintf(fgseaRes$NES[2], fmt = '%#.2f'), "\npadj: ", 
                     sprintf(fgseaRes$padj[2], fmt = '%#.2e')))
    
    print(p.down)
    dev.off()
    
    message('===========================\n')
  }
}

  