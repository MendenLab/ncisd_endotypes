rm(list = ls(all=TRUE))

# library(lme4)
library(variancePartition, lib.loc='/Users/christina.hillig/R_packages')  # BiocManager::install("variancePartition", lib='/Users/christina.hillig/R_packages')
library(edgeR)
library(limma)
library(DESeq2)
library('xlsx')


# To estimate the variance in lesional gene expression explained by differences in 
# diagnosis, gender, age and other clinical attributes, 
# the RNAseq count data was first transformed using voom transformation from R package Limma and then modelled with a linear mixed effect model using bioconductor package variancePartition.


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

# 1. Data preparation
obs.name <- "healthysamp_diag_v2"
# 1.1 Load data
dds <- def.load_data(input.dir="/Users/christina.hillig/Documents/Projects/Eyerich_AG_projects/BRAIN__Peter_Seiringer")
# Load suptypes
input.dir.subtypes <- '/Users/christina.hillig/PycharmProjects/Eyerich_ERCGrant/Molecular_subtypes/output/Leiden/compare_resolutions/2023-02-21'
df <- xlsx::read.xlsx(file.path(input.dir.subtypes, 'Optimalres_0.5__OptimalGPTK_7.xlsx'), sheetIndex = 1)
subtypes.list.res0.5 <- def.prepare_subtypes(df.file = df)
subtypes.list.res0.5 = as.integer(subtypes.list.res0.5)
# drop NAs
subtypes.list.res0.5 <- subtypes.list.res0.5[!is.na(subtypes.list.res0.5)]
df <- xlsx::read.xlsx(file.path(input.dir.subtypes, 'Optimalres_0.8__OptimalGPTK_7.xlsx'), sheetIndex = 1)
subtypes.list.res0.8 <- def.prepare_subtypes(df.file = df)
subtypes.list.res0.8 = as.integer(subtypes.list.res0.8)
# drop NAs
subtypes.list.res0.8 <- subtypes.list.res0.8[!is.na(subtypes.list.res0.8)]
df <- xlsx::read.xlsx(file.path(input.dir.subtypes, 'Optimalres_0.9__OptimalGPTK_7.xlsx'), sheetIndex = 1)
subtypes.list.res0.9 <- def.prepare_subtypes(df.file = df)
subtypes.list.res0.9 = as.integer(subtypes.list.res0.9)
# drop NAs
subtypes.list.res0.9 <- subtypes.list.res0.9[!is.na(subtypes.list.res0.9)]

# Remove samples which are not in the only lesion object ..
mucids.L <- def.prepare_subtypes(df.file=df, obs='MUC.IDs') 
mucids.NL <- dds[,(dds$sampleType == 'H')]$sampleID
mucids <- c(mucids.L, mucids.NL)
dds <- dds[ , dds$sampleID %in% mucids]
# # Find matching MUC IDs
ind.dds <- match(mucids.L, colnames(dds[ , dds$sampleType == 'D']))

# Add subtypes
dds$`Molecular subtype res0.5` <- 'NL'
dds$`Molecular subtype res0.8` <- 'NL'
dds$`Molecular subtype res0.9` <- 'NL'
dds$`Molecular subtype res0.5`[dds$sampleType == 'D'][ind.dds] <- subtypes.list.res0.5
dds$`Molecular subtype res0.8`[dds$sampleType == 'D'][ind.dds] <- subtypes.list.res0.8
dds$`Molecular subtype res0.9`[dds$sampleType == 'D'][ind.dds] <- subtypes.list.res0.9


# Create edgeR object
y.object <- def.create.dgelistobject(
  counts=counts(dds), metadata=colData(dds), genes.infos=rowData(dds), 
  group.infos=dds[[obs.name]]) 
# Perform TMM normalization
y.object <- calcNormFactors(y.object)

# 1.2 Create design matrix
design.matrix <- model.matrix(~ 1, data = y.object$samples)

# 1.3 Filter counts
y.object = def.prepare_dgelistobject(y.obj=y.object)

# 1.4 Transform counts
v.object <- limma::voom(y.object, design.matrix, plot=TRUE)

# 2. variancePartition
# 2.1 Formula indicating which meta-data variables to consider 
formula.design = formula(~ (1|diag) + (1|batchID) + age + (1|Sex.x) + (1|sampleType) + 
                           (1|Molecular.subtype.res0.5) + 
                           (1|Molecular.subtype.res0.8) + (1|Molecular.subtype.res0.9))

# 2.2 MetaData
info <- data.frame(diag = y.object$samples$diag, batchID = y.object$samples$batchID,
                   age = y.object$samples$age, Sex.x= y.object$samples$Sex.x,
                   Pattern = y.object$samples$Pattern, sampleType= y.object$samples$sampleType,
                   Molecular.subtype.res0.5 = y.object$samples$Molecular.subtype.res0.5, 
                   Molecular.subtype.res0.8 = y.object$samples$Molecular.subtype.res0.8, 
                   Molecular.subtype.res0.9 = y.object$samples$Molecular.subtype.res0.9)

# 2.2 Run variancePartition
varPart <- variancePartition::fitExtractVarPartModel(
  exprObj=v.object, formula=formula.design, data=info )

# TODO @Guelce: please save the varPart object so that I can use it for further analysis :)

# vp <- sortCols( varPart )
# plotPercentBars( vp[1:10,] )
# plotVarPart( vp )


# 3. Removing batch effects prior to variancePartition to assess the
# influence of biological variables better
# subtract out effect of Batch with linear mixed model
# Too memory expesive 
# modelFit <- variancePartition::fitVarPartModel(
#   exprObj=v.object, formula=~ (1|batchID) + age + (1|Sex.x), data=info ) 
# residMatrix <- variancePartition::residuals( modelFit )
# Mempry efficient: extract residuals directly without storing intermediate results
residList <- fitVarPartModel( exprObj=v.object, formula=~ (1|batchID) + age + (1|Sex.x), data=info )
# convert list to matrix
residMatrix = do.call(rbind, residList)

# fit model on residuals
form.v2 <- formula(~ (1|diag) + (1|Molecular.subtype.res0.5) + (1|sampleType) +
                     (1|Molecular.subtype.res0.8) + (1|Molecular.subtype.res0.9))
varPartResid <- variancePartition::fitExtractVarPartModel(
  exprObj=residMatrix, formula=form.v2, data=info )

# TODO @Guelce: please save also the varPartResid object so that I can use it for further analysis :)


