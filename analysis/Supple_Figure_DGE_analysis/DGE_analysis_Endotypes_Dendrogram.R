# Script to perform DGE analysis between branches of a node
rm(list = ls(all=TRUE))

library(edgeR)
library(limma)
library(DESeq2)
library('xlsx')
library(openxlsx)
library(data.table)
library("dendextend")

root.path <- '/Volumes/CH__data/Projects/Eyerich_AG_projects/BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer'
script.path <- file.path(root.path, 'analysis', 'Molecular_subtypes', 'r_scripts', 'DGE_analysis')
source(file.path(script.path, 'utils.R'))
source(file.path(script.path, 'Volcanoplot.R'))
source(file.path(script.path, 'perform_DGE_analysis.R'))


#' @description Load .rds object containing lesion and non-lesion samples including metaData
#' @param input.dir Path to count matrix and meta.data
def.load_data <- function(input.dir, name.dds.obj='dds_highQual_Sexandbatchcorrected_v04.rds') 
{
  
  # Load raw expression data but filtered for protein transcribing genes and TPM > 1
  # Design: ~ SampleType + Sex.x + batch
  dds <- readRDS(file.path(input.dir, 'count_matrix', 'data_freeze', 'preprocessed', 'L',
                           name.dds.obj))
  # dim(counts(dds)): 18341   632
  
  # Load metaData
  meta.data <- xlsx::read.xlsx2(
    file.path(input.dir, 'clinical_data', '20210720_patient_meta_data_v04__CH.xlsx'), sheetIndex = 1)
  
  ind.samples <- match(dds$sampleID, meta.data$Helmholtz.identifyer)
  meta.data <- meta.data[ind.samples, ]
  row.names(meta.data) <- meta.data$Helmholtz.identifyer
  
  testit::assert(rownames(colData(dds)) == row.names(meta.data))
  
  
  df.temp <- colData(dds)[c('Year_sequencing', 'RIN_D', 'DV200_D', 'healthysamp_diag', 
                            'healthysamp_pattern', 'batch_merged', 'sizeFactor', 'replaceable')] 
  # Overwrite batchID as it did not contain all values ..
  dds$batchID <- meta.data$batchID
  dds$batchID <- as.factor(dds$batchID)
  
  # # overwrite dds$healthysamp_diag as it contains NA values ...
  # dds <- def.rename_nonlesion_samples_to_paired_lesion_name(
  #   dds=dds, ref.obs.name='diag', obs.name='healthysamp_diag_v2', group.name = 'diag') 
  
  # add ENSEMBL to rowData
  rowData(dds)$ENSEMBL <- rowData(dds)$ensembl_id
  rowData(dds)$entrezid <- rowData(dds)$entrezgene_id
  rowData(dds)$hgnc_symbol <- make.unique(rowData(dds)$Gene_name)
  rowData(dds)$description <- rowData(dds)$entrezgene_description
  rowData(dds)$gene_name <- rowData(dds)$entrezgene_description
  
  # Rename rows to gene id (here: entrezgene_accession)
  row.names(dds) <- make.unique(rowData(dds)$Gene_name)
  
  return(dds)
  
}


def.create.dgelistobject <- function(counts, metadata, genes.infos, group.infos) 
{
  y <- edgeR::DGEList(counts = counts, samples = metadata, genes = genes.infos, 
                      group = group.infos, remove.zeros = TRUE)
  
  
  return(y)
}


def.prepare_dgelistobject <- function(y.obj) 
{
  # filter count matrix
  message('Keep only counts >= 0')
  
  # min.total.count = 2 so that gene was at least measured in two samples
  keep <- edgeR::filterByExpr(y.obj, min.count = 1, min.total.count = 2)
  y.obj <- y.obj[keep, , keep.lib.sizes=FALSE]
  
  return(y.obj)
}


# 0. Parameters 
obs.name <- "healthysamp_clusterlabels"
opt.method = 'cosine'
opt.linkage_method = 'ward.D2'

output.dir <- file.path(root.path, 'analysis', 'Molecular_subtypes', 'output')
output.dir <- file.path(output.dir, 'DGE_analysis', 'Endotypes_Dendrogram', Sys.Date())
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

# 1. Load data
df.grouped.counts <- read.csv(file.path(
  root.path, 'analysis', 'Molecular_subtypes', 'output',
  'Dendrogram', '2023-09-13', 'Normed_Dataframe_mean_gene_selection.csv'))
# '2023-05-31'
# Reorder rows by cluster names
sorting <- order(as.integer(gsub("E", "\\1", df.grouped.counts$cluster)))
df.grouped.counts <- df.grouped.counts[ sorting, ]
cluster.names <- df.grouped.counts$cluster
row.names(df.grouped.counts) <- seq(1, dim(df.grouped.counts)[1])
df.grouped.counts$cluster <- NULL

dds <- def.load_data(
  input.dir=file.path(root.path, 'raw_data'), 
  name.dds.obj='dds_highQual_Sexandbatchcorrected__Endotypes_renamed_v01.rds')
# subset to lesion samples
dds <- dds[, dds$sampleType != 'H']

# Create edgeR object
y.object <- def.create.dgelistobject(
  counts=counts(dds), metadata=colData(dds), genes.infos=rowData(dds), 
  group.infos=dds$Molecular_subtype) 


# 2. Create dendrogram
Matrix <- as.matrix(df.grouped.counts)
res.dist <- Matrix / sqrt(rowSums(Matrix * Matrix))
res.dist <- res.dist %*% t(res.dist)
res.dist <- as.dist(1 - res.dist)
Z <- hclust(d=res.dist, method=opt.linkage_method)
tree <- ape::as.phylo(Z)

# 3. Create DGE analysis pairs
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



# 4. Run DGE analysis one vs one
message('DGE analysis of L skin biopsies within branches of Dendrogram')
hash.dict <- hash::hash()

for (row_number in seq(1, dim(df.pairs)[1])) 
{
  df.comb.temp <- df.pairs[as.integer(row_number), ]
  val.1 <- as.integer(unlist(strsplit(df.comb.temp[1],", "))) * -1
  val.2 <- as.integer(unlist(strsplit(df.comb.temp[2],", "))) * -1
  
  val.1 <- paste("E", val.1, sep='')
  val.2 <- paste("E", val.2, sep='')
  
  
  message(paste0(paste(val.1, collapse = '_'), ' vs ',
                 paste(val.2, collapse = '_')))
  
  y_subset <- y.object[, y.object$samples$Molecular_subtype %in% val.1 | 
                         y.object$samples$Molecular_subtype %in% val.2]
  # Add new column for comparison
  y_subset$samples$comparison <- 0
  # Select Rows by list of column Values
  y_subset$samples$comparison[y_subset$samples$Molecular_subtype %in% val.1] <- 1
  y_subset$samples$comparison[y_subset$samples$Molecular_subtype %in% val.2] <- 2
  # Convert to factor
  y_subset$samples$comparison <- as.factor(y_subset$samples$comparison)
  
  # relevel so that reference (2) is always first 
  y_subset$samples$comparison <- stats::relevel(
    y_subset$samples$comparison , ref = 2)
  
  y_subset <- def.prepare_dgelistobject(y.obj=y_subset) 
  
  
  tryCatch(
    {
      design.name = 'Age_Sex_batchID_subtypes'
      comparison.name <- paste(paste(val.1, collapse = '_'), ' vs ',
                               paste(val.2, collapse = '_'), 
                               design.name, sep = "__")
      
      design.matrix <- model.matrix(
        ~age + Sex.x + batchID + comparison, data = y_subset$samples)
      
      # 4. Run DGE analysis using EdgeR
      dge.result_dge.object <- def.run_general_edger(
        y.obj=y_subset, design.matrix = design.matrix) 
      
      
    }, error=function(e)
    {
      design.name = 'Age_Sex_subtypes'
      comparison.name <- paste(paste0(paste(val.1, collapse = '_'), ' vs ',
                                      paste(val.2, collapse = '_')), 
                               design.name, sep = "__")
      
      message('Switch to design matrix without correcting for batchID')
      design.matrix <- model.matrix(
        ~age + Sex.x + comparison, data = y_subset$samples)
      
      # 4. Run DGE analysis using EdgeR
      dge.result_dge.object <- def.run_general_edger(
        y.obj=y_subset, design.matrix = design.matrix) 
    })
  
  hash.dict[[paste0(paste(val.1, collapse = '_'), ' vs ',
                    paste(val.2, collapse = '_'))]] <- dge.result_dge.object[[1]]
  
  message("Save to csv and xlsx file")
  write.csv(dge.result_dge.object[[1]],
            file.path(output.dir, paste0('DEGs_', comparison.name, '.csv')))
  openxlsx::write.xlsx(dge.result_dge.object[[1]],
             file.path(output.dir, paste0('DEGs_', comparison.name, '.xlsx')))
  
  
  # 5. Plot Volcano with marker genes
  message("Volcanoplot")
  pdf(file = file.path(output.dir, paste0(comparison.name, '_Volcanoplot', '.pdf')), 
      width = 6, height = 6, family = "Helvetica")
  print(
    def.get.volcanoplot(df.data=dge.result_dge.object[[1]], log2fc.cut=1, 
                        add.labels=NULL, add.annotations=TRUE) 
  )
  dev.off()
}

