rm(list = ls(all=TRUE))

library(edgeR)
library(limma)
library(DESeq2)
library('xlsx')
library(data.table)

source('/Users/christina.hillig/R_studio/K_ERC_Grant/scripts/DGE_analysis/Volcanoplot.R')
source('/Users/christina.hillig/R_studio/K_ERC_Grant/scripts/DGE_analysis/perform_DGE_analysis.R')


#' @description Load .rds object containing lesion and non-lesion samples including metaData
#' @param input.dir Path to count matrix and meta.data
def.load_data <- function(input.dir) 
{
  
  # Load raw expression data but filtered for protein transcribing genes and TPM > 1
  # Design: ~ SampleType + Sex.x + batch
  dds <- readRDS(file.path(input.dir, 'count_matrix', 'data_freeze', 'preprocessed', 'L_and_NL',
                           'dds_highQual_Sexandbatchcorrected_v04.rds'))
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
  
  # overwrite dds$healthysamp_diag as it contains NA values ...
  dds <- def.rename_nonlesion_samples_to_paired_lesion_name(
    dds=dds, ref.obs.name='diag', obs.name='healthysamp_diag_v2') 
  
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


def.rename_nonlesion_samples_to_paired_lesion_name <- function(dds, ref.obs.name, obs.name) 
{
  dds[[ref.obs.name]] <- as.factor(dds[[ref.obs.name]])
  
  levels.ref.obs <- levels(dds[[ref.obs.name]])
  # remove non-lesional to set it on first place
  levels.ref.obs <- levels.ref.obs[levels.ref.obs != "non-lesional"]
  # Reorder levels such that non-lesional level comes first
  dds[[ref.obs.name]] <- factor(dds[[ref.obs.name]], levels=c('non-lesional', levels.ref.obs))
  
  print(levels(dds[[ref.obs.name]]))
  
  # overwrite dds$healthysamp_diag as it contains NA values ...
  dds[[obs.name]] <- 'Unknown'
  # 40 only non-lesion samples (after filtering and everything)
  for (observable in levels(dds[[ref.obs.name]]))
  {
    # print(observable)
    mask.observable <- dds[[ref.obs.name]] == observable
    # Find partner via Patient ID == Pseudo.ID
    patient.id <- dds$PatientID[mask.observable]
    # multiple matching using which and %in%, as match only return first encounter 
    mask.patients <- which(dds$PatientID %in% patient.id)
    dds[[obs.name]][mask.patients] <- observable
  }
  dds[[obs.name]] <- as.factor(dds[[obs.name]])
  
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

def.prepare_subtypes <- function(df.file, obs='Molecular.subtypes') 
{
  # Save subtypes for seed 0
  subtypes.list <- df.file[1, obs]
  subtypes.list <- (gsub("\\[|\\]|[\r\n]|'|,", "", subtypes.list))
  subtypes.list <-  strsplit(subtypes.list, " ")[[1]]
  
  return(subtypes.list)
}

output.dir <- '/Users/christina.hillig/R_studio/K_ERC_Grant/output'
output.dir <- file.path(output.dir, 'DGE_analysis', 'Molecular_subtypes_res0.9', Sys.Date())
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

# 1. Data preparation
obs.name <- "healthysamp_diag_v2"
# 1.1 Load data
dds <- def.load_data(input.dir="/Users/christina.hillig/Documents/Projects/Eyerich_AG_projects/BRAIN__Peter_Seiringer")
# Remove Non Lesion samples
dds <- dds[ , dds$sampleType == 'D']
# Load suptypes
input.dir.subtypes <- '/Users/christina.hillig/PycharmProjects/Eyerich_ERCGrant/Molecular_subtypes/output/Leiden/compare_resolutions/2023-02-21'
df <- xlsx::read.xlsx(file.path(input.dir.subtypes, 'Optimalres_0.9__OptimalGPTK_7.xlsx'), sheetIndex = 1)
subtypes.list.res0.9 <- def.prepare_subtypes(df.file = df)
subtypes.list.res0.9 = as.integer(subtypes.list.res0.9)
# drop NAs
subtypes.list.res0.9 <- subtypes.list.res0.9[!is.na(subtypes.list.res0.9)]

# Remove samples which are not in the only lesion object ..
mucids.L <- def.prepare_subtypes(df.file=df, obs='MUC.IDs') 
dds <- dds[ , dds$sampleID %in% mucids.L]
# # Find matching MUC IDs
ind.dds <- match(mucids.L, colnames(dds[ , dds$sampleType == 'D']))

# Add subtypes
dds$`Molecular_subtype` <- 'L'
dds$`Molecular_subtype`[dds$sampleType == 'D'][ind.dds] <- subtypes.list.res0.9
dds$`Molecular_subtype` <- as.factor(dds$`Molecular_subtype`)

# Create edgeR object
y.object <- def.create.dgelistobject(
  counts=counts(dds), metadata=colData(dds), genes.infos=rowData(dds), 
  group.infos=dds[[obs.name]]) 


# 2. Run DGE analysis one vs one
# 3. Perform DGE analysis between each pattern 
message('DGE analysis of L skin biopsies between Endotypes')
hash.dict <- hash::hash()
table(y.object$samples$Molecular_subtype, y.object$samples$Molecular_subtype)

df.combinations <- data.table::CJ(levels(y.object$samples$Molecular_subtype), 
                                  levels(y.object$samples$Molecular_subtype))  
df.combinations <- df.combinations[!duplicated(t(apply(df.combinations, 1, sort))), ]

for (row_number in row.names(df.combinations)) 
{
  df.comb.temp <- df.combinations[as.integer(row_number), ]
  if (df.comb.temp$V1 != df.comb.temp$V2) 
  {
    message(paste0(df.comb.temp$V1, '_vs_', df.comb.temp$V2))
    y_subset <- y.object[, y.object$samples$Molecular_subtype == df.comb.temp$V1 | 
                      y.object$samples$Molecular_subtype == df.comb.temp$V2]
    y_subset$samples$Molecular_subtype <- droplevels(y_subset$samples$Molecular_subtype)
    
    if (dim(y_subset)[2] > 3 & length(unique(y_subset$samples$Molecular_subtype)) >=2 & 
        min(table(y_subset$samples$Molecular_subtype)) >=3) 
    {
      # relevel so that reference is always first 
      y_subset$samples$Molecular_subtype <- stats::relevel(
        y_subset$samples$Molecular_subtype , ref = df.comb.temp$V2)
      # Create contrast - Alphabetical order!!
      contrast <- c('condition', levels(y_subset$samples$Molecular_subtype)[1], 
                    levels(y_subset$samples$Molecular_subtype)[2])
      
      y_subset <- def.prepare_dgelistobject(y.obj=y_subset) 
      
      # Important: Create design matrix after releveling
      # TODO add try catch because
      # Error in glmFit.default(sely, design, offset = seloffset, dispersion = 0.05,  : 
      # Design matrix not of full rank.  The following coefficients not estimable: batchID7
      # full model matrix is less than full rank if I include batch ..
      tryCatch(
        {
          design.name = 'Age_Sex_batchID_Molecular_subtype'
          comparison.name <- paste(paste0(df.comb.temp$V1, '_vs_', df.comb.temp$V2), design.name, sep = "__")
          
          design.matrix <- model.matrix(
            ~age + Sex.x + batchID + Molecular_subtype, data = y_subset$samples)
          
          # 4. Run DGE analysis using EdgeR
          dge.result_dge.object <- def.run_general_edger(
            y.obj=y_subset, contrast=contrast, design.matrix = design.matrix) 
          

        }, error=function(e)
        {
          design.name = 'Age_Sex_Molecular_subtype'
          comparison.name <- paste(paste0(df.comb.temp$V1, '_vs_', df.comb.temp$V2), design.name, sep = "__")
          
          message('Switch to design matrix without correcting for batchID')
          design.matrix <- model.matrix(
            ~age + Sex.x + Molecular_subtype, data = y_subset$samples)
          
          # 4. Run DGE analysis using EdgeR
          dge.result_dge.object <- def.run_general_edger(
            y.obj=y_subset, contrast=contrast, design.matrix = design.matrix) 
        })
      
      hash.dict[[
        paste0(df.comb.temp$V1, '_vs_', df.comb.temp$V2)]] <- dge.result_dge.object[[1]]
      
      write.csv(dge.result_dge.object[[1]],
                file.path(output.dir, paste0('DEGs_', comparison.name, '.csv')))
      write.xlsx(dge.result_dge.object[[1]],
                file.path(output.dir, paste0('DEGs_', comparison.name, '.xlsx')))
      
      
      # 5. Plot Volcano with marker genes
      pdf(file = file.path(output.dir, paste0(comparison.name, '_Volcanoplot', '.pdf')), 
          width = 6, height = 6)
      print(
        def.get.volcanoplot(df.data=dge.result_dge.object[[1]], log2fc.cut=1, add.labels=NULL) 
      )
      dev.off()
      
    }
  }
}
