rm(list = ls(all=TRUE))

# library(lme4)
# library(variancePartition, lib.loc='/Users/christina.hillig/R_packages')  # BiocManager::install("variancePartition", lib='/Users/christina.hillig/R_packages')
library(edgeR)
library(limma)
library(DESeq2)
library('xlsx')
library(data.table)

root.path <- '/Volumes/CH__data/Projects/Eyerich_AG_projects/BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer'
script.path <- file.path(root.path, 'analysis', 'Molecular_subtypes', 'r_scripts', 'DGE_analysis')
source(file.path(script.path, 'Volcanoplot.R'))
source(file.path(script.path, 'perform_DGE_analysis.R'))


#' @description Load .rds object containing lesion and non-lesion samples including metaData
#' @param input.dir Path to count matrix and meta.data
def.load_data <- function(input.dir, name.dds.obj='dds_highQual_Sexandbatchcorrected_v04.rds') 
{
  
  # Load raw expression data but filtered for protein transcribing genes and TPM > 1
  # Design: ~ SampleType + Sex.x + batch
  dds <- readRDS(file.path(input.dir, 'count_matrix', 'data_freeze', 'preprocessed', 'L_and_NL',
                           name.dds.obj))
  # dim(counts(dds)): 18341   632
  
  # Load metaData
  meta.data <- xlsx::read.xlsx2(
    file.path(input.dir, 'clinical_data', '20210720_patient_meta_data_v04__CH.xlsx'), sheetIndex = 1)
  
  ind.samples <- match(dds$sampleID, meta.data$Helmholtz.identifyer)
  meta.data <- meta.data[ind.samples, ]
  row.names(meta.data) <- meta.data$Helmholtz.identifyer
  
  testit::assert(rownames(SummarizedExperiment::colData(dds)) == row.names(meta.data))
  
  
  df.temp <- SummarizedExperiment::colData(dds)[
    c('Year_sequencing', 'RIN_D', 'DV200_D', 'healthysamp_diag', 
      'healthysamp_pattern', 'batch_merged', 'sizeFactor', 'replaceable')] 
  # Overwrite batchID as it did not contain all values ..
  dds$batchID <- meta.data$batchID
  dds$batchID <- as.factor(dds$batchID)
  
  # overwrite dds$healthysamp_diag as it contains NA values ...
  dds <- def.rename_nonlesion_samples_to_paired_lesion_name(
    dds=dds, ref.obs.name='diag', obs.name='healthysamp_diag_v2') 
  
  # add ENSEMBL to rowData
  SummarizedExperiment::rowData(dds)$ENSEMBL <- SummarizedExperiment::rowData(dds)$ensembl_id
  SummarizedExperiment::rowData(dds)$entrezid <- SummarizedExperiment::rowData(dds)$entrezgene_id
  SummarizedExperiment::rowData(dds)$hgnc_symbol <- make.unique(SummarizedExperiment::rowData(dds)$Gene_name)
  SummarizedExperiment::rowData(dds)$description <- SummarizedExperiment::rowData(dds)$entrezgene_description
  SummarizedExperiment::rowData(dds)$gene_name <- SummarizedExperiment::rowData(dds)$entrezgene_description
  
  # Rename rows to gene id (here: entrezgene_accession)
  row.names(dds) <- make.unique(SummarizedExperiment::rowData(dds)$Gene_name)
  
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

output.dir <- file.path(root.path, 'analysis', 'Molecular_subtypes', 'output')
output.dir <- file.path(output.dir, 'DGE_analysis', 'Molecular_subtypes_res0.9_1vsRest', Sys.Date())
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

# 1. Data preparation
obs.name <- "healthysamp_diag_v2"
# 1.1 Load data
dds <- def.load_data(
  input.dir=file.path(root.path, 'raw_data'), 
  name.dds.obj='dds_highQual_Sexandbatchcorrected__Endotypes_renamed_v01.rds')
# 1.2 subset to lesion samples
dds <- dds[, dds$sampleType != 'H']


# Create edgeR object
# Create edgeR object
y.object <- def.create.dgelistobject(
  counts=counts(dds), metadata=SummarizedExperiment::colData(dds), 
  genes.infos=SummarizedExperiment::rowData(dds), 
  group.infos=dds$Molecular_subtype) 


# 2. Run DGE analysis cluster vs rest
# 3. SPerform DGE analysis between each pattern 
message('DGE analysis of L skin biopsies between Endotypes')
hash.dict <- hash::hash()
table(y.object$samples$Molecular_subtype, y.object$samples$Molecular_subtype)

df.combinations <- data.table::CJ(levels(y.object$samples$Molecular_subtype), 
                                  levels(y.object$samples$Molecular_subtype))  
df.combinations <- df.combinations[!duplicated(t(apply(df.combinations, 1, sort))), ]

# Drop unused factors
y.object$samples$Molecular_subtype <- factor(y.object$samples$Molecular_subtype)

subtypes.list <- levels(y.object$samples$Molecular_subtype)

for (subtype.name in subtypes.list) 
{
  # TODO 
  y.object$samples$comparision <- subtype.name
  mask.tmp <- y.object$samples$Molecular_subtype != subtype.name
  y.object$samples$comparision[mask.tmp] <- "Rest"

  message(paste0(subtype.name, '_vs_', "Rest"))
  design.name = 'Age_Sex_batchID_Molecular_subtype'
  comparison.name <- paste(paste0(subtype.name, '_vs_', "Rest"), 
                           design.name, sep = "__")
  y.object$samples$comparision <- as.factor( y.object$samples$comparision)
  
  if (dim(y.object)[2] > 3 & length(unique(y.object$samples$comparision)) >=2 & 
      min(table(y.object$samples$comparision)) >=3) 
  {
    # relevel so that reference is always first 
    y.object$samples$comparision <- stats::relevel(
      y.object$samples$comparision , ref = "Rest")
    # Create contrast - Alphabetical order!!
    contrast <- c('condition', levels(y.object$samples$comparision)[1], 
                  levels(y.object$samples$comparision)[2])
    # Important: Create design matrix after releveling
    design.matrix <- model.matrix(~age + Sex.x + batchID + comparision, 
                                  data = y.object$samples)
    
    y_subset <- def.prepare_dgelistobject(y.obj=y.object) 
    
    # 4. Run DGE analysis using EdgeR
    dge.result_dge.object <- def.run_general_edger(
      y.obj=y_subset, design.matrix = design.matrix) 
    
    hash.dict[[paste0(subtype.name, '_vs_', "Rest")]] <- dge.result_dge.object[[1]]
    
    write.csv(dge.result_dge.object[[1]],
              file.path(output.dir, paste0('DEGs_', comparison.name, '.csv')))
    write.xlsx(dge.result_dge.object[[1]],
               file.path(output.dir, paste0('DEGs_', comparison.name, '.xlsx')))
    
    
    # 5. Plot Volcano with marker genes
    pdf(file = file.path(output.dir, paste0(comparison.name, '_Volcanoplot', '.pdf')), 
        width = 6, height = 6, family = "Helvetica")
    print(
      def.get.volcanoplot(df.data=dge.result_dge.object[[1]], log2fc.cut=1, 
                          add.labels=NULL) 
    )
    dev.off()
    
  }
}
