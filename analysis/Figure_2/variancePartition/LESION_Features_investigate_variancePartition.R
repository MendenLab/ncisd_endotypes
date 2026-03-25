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
                            'healthysamp_pattern', 'batch_merged', 'sizeFactor',
                            'replaceable')] 
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


######################################################################
#                                Start                               #
######################################################################

# 1. Data preparation
obs.name <- "healthysamp_diag_v2"
# 1.1 Load data
dds <- def.load_data(input.dir="/Users/christina.hillig/Documents/Projects/Eyerich_AG_projects/BRAIN__Peter_Seiringer")
# Remove Non Lesion samples
dds <- dds[ , dds$sampleType == 'D']


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
dds <- dds[ , dds$sampleID %in% mucids.L]
# # Find matching MUC IDs
ind.dds <- match(mucids.L, colnames(dds[ , dds$sampleType == 'D']))

# Add subtypes to DESeq2 object
dds$`Molecular subtype res0.5` <- 'NL'
dds$`Molecular subtype res0.8` <- 'NL'
dds$`Molecular subtype res0.9` <- 'NL'
dds$`Molecular subtype res0.5`[dds$sampleType == 'D'][ind.dds] <- subtypes.list.res0.5
dds$`Molecular subtype res0.8`[dds$sampleType == 'D'][ind.dds] <- subtypes.list.res0.8
dds$`Molecular subtype res0.9`[dds$sampleType == 'D'][ind.dds] <- subtypes.list.res0.9


# group clinical attributes by Age, Sex, Morphology, Severity, Therapy, Histology, Comorbidity, and History 
cattr_age <- 'age.x'
cattr_sex <- 'Sex.x'
cattr_morphology <- c('Morph_erythema','Morph_Elev', 'Morph_papu', 'Morph_vesicles',
                      'Morph_pustules', 'Morph_scaling', 'Morph_dryness')
cattr_severity <- c('Diag_PGA.x', 'Diag_PASI', 'Diag_SCORAD', 'Diag_DLQI')
cattr_therapy <- c('Therapeutic.target_IL.17', 'Therapeutic.target_IL.23', 
                   'Therapeutic.target_MTX', 'Therapeutic.target_TNF', 
                   'The_dor', 'The_sor', 'The_leuco', 'The_Hb', 'The_lympho', 'The_weeks',
                   'The_granulo', 'The_eosino', 'The_crea', 'The_GPT')
cattr_laboratory <- c('Lab_csa', 'Lab_specificIgE', 'Lab_leuco', 'Lab_Hb', 'Lab_lympho',
                      'Lab_granulo', 'Lab_eosino', 'Lab_crea', 'Lab_GPT', 'Lab_totalIgE')
cattr_histology <- c('Hist_SE', 'Hist_Para_quali', 'Hist_Aka_quali', 'Hist_Granu',
                     'Hist_Serum',
                     'Hist_Dist_Lympho', 'Hist_Microabscess', 'Hist_Bacteria',
                     'Hist_Neutro_quali', 'Hist_Keratoses', 'Hist_ID_distribution',
                     'Hist_ID_quali', 'Hist_Nr_Keratoses', 'Hist_Hyper_quant',
                     'Hist_Para_quant',
                     'Hist_Ortho_quant', 'Hist_Aka_quant', 'Hist_Spongiosis', 'Hist_Cap',
                     'Hist_Lympho', 'Hist_Exocytosis', 'Hist_Neutro_quant', 
                     'Hist_Eosino', 'Hist_ID_quant', 'Hist_LCV', 'Hist_Mucin')
cattr_comorbidity <- c('Com_arths', 'Com_ah', 'Com_iaa', 'Com_rca', 'Com_dm', 'Com_ast', 
                       'Com_hepa', 'Com_renal', 'Com_tonsil', 'Com_ibd', 'Com_nail', 
                       'Com_muc', 'Com_scalp', 'Com_BMI')
cattr_history <- c('Hty_fh', 'Hty_exanth_1.0',
                   'Hty_exanth_2.0', 'Hty_self', 'Hty_photos', 'Hty_progr', 'Hty_smoker',
                   'Hty_sport', 'Hty_alc', 'Hty_onset', 'Hty_pruritus', 'Hty_chronic')
## add clinical attributes
df.ca <- read.csv('/Users/christina.hillig/PycharmProjects/Eyerich_ERCGrant/Molecular_subtypes/output/Figure_1B_Heatmap/Clinical_attributes.csv')
df.ca$sampleID <- df.ca$Helmholtz_identifyer
row.names(df.ca) <- df.ca$Helmholtz_identifyer
# Merge metaData to dds object
df.tmp <- as.data.frame(merge(x = colData(dds), y = df.ca, by = "sampleID"))
row.names(df.tmp) <- df.tmp$sampleID
dds = DESeq2::DESeqDataSetFromMatrix(
  countData=counts(dds), colData=df.tmp, design=~age.x + batchID + Sex.x + diag)


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


######################################################################
#         variancePartition without removing batch effects           #
######################################################################
columns.to.investigate <- c('diag', 'Pattern', 'batchID', 'sampleType',
                            'Molecular.subtype.res0.5', 'Molecular.subtype.res0.8',
                            'Molecular.subtype.res0.9',
                            cattr_age, cattr_sex, cattr_morphology, cattr_severity,
                            cattr_therapy, cattr_laboratory, cattr_histology,
                            cattr_comorbidity, cattr_history)
factors.to.investigate <- c(
  'Morph_vesicles', 'Diag_PGA.x', 'Diag_SCORAD', 'Therapeutic.target_IL.17', 
  'Therapeutic.target_IL.23', 'Therapeutic.target_MTX', 'Therapeutic.target_TNF', 
  'Lab_csa', 'Lab_specificIgE', 'Hist_Bacteria', 'Hist_Keratoses', 'Hist_ID_quant',
  'Hist_LCV', 'Com_arths', 'Com_ah', 'Com_iaa', 'Com_rca', 'Com_dm', 'Com_ast',
  'Com_hepa', 'Com_renal', 'Com_tonsil', 'Com_ibd', 'Com_nail', 'Com_muc', 
  'Com_scalp', 'Hty_fh', 'Hty_exanth_1.0', 'Hty_exanth_2.0', 'Hty_self','Hty_photos',
  'Hty_progr', 'Hty_smoker', 'Hty_onset' )
for (attr.name in columns.to.investigate) 
{
  y.object$samples[[attr.name]] <- as.factor(y.object$samples[[attr.name]] )
}

# 2. variancePartition
# 2.1 MetaData
info <- y.object$samples %>% select(columns.to.investigate)

# 2.3 Run variancePartition 
# 2.3.1 res 0.5 
formula.design = formula(
  ~ (1|Molecular.subtype.res0.5) + (1|diag) + (1|batchID) +
    (1|Pattern)  +  age.x + (1|Sex.x) +
    Morph_erythema + Morph_Elev + Morph_papu + (1|Morph_vesicles) + Morph_pustules +
    Morph_scaling + Morph_dryness + (1|Diag_PGA.x) + Diag_PASI + (1|Diag_SCORAD) + 
    Diag_DLQI + (1|Therapeutic.target_IL.17) + (1|Therapeutic.target_IL.23) + 
    (1|Therapeutic.target_MTX) + (1|Therapeutic.target_TNF) + The_dor + The_sor + 
    The_leuco + The_Hb + The_lympho + The_weeks + The_granulo + The_eosino + The_crea +
    The_GPT + (1|Lab_csa) + (1|Lab_specificIgE) + Lab_leuco + Lab_Hb + Lab_lympho +
    Lab_granulo + Lab_eosino + Lab_crea + Lab_GPT + Lab_totalIgE +  Hist_SE + 
    Hist_Para_quali + Hist_Aka_quali + Hist_Granu + Hist_Serum + Hist_Dist_Lympho +
    Hist_Microabscess + (1|Hist_Bacteria) + Hist_Neutro_quali + (1|Hist_Keratoses) +
    Hist_ID_distribution + Hist_ID_quali + Hist_Nr_Keratoses + Hist_Hyper_quant +
    Hist_Para_quant + Hist_Ortho_quant + Hist_Aka_quant + Hist_Spongiosis + 
    Hist_Cap + Hist_Lympho + Hist_Exocytosis + Hist_Neutro_quant + Hist_Eosino +
    (1|Hist_ID_quant) + (1|Hist_LCV) + Hist_Mucin + (1|Com_arths) + (1|Com_ah) + 
    (1|Com_iaa) + (1|Com_rca) + (1|Com_dm) + (1|Com_ast) + (1|Com_hepa) + (1|Com_renal) + 
    (1|Com_tonsil) + (1|Com_ibd) + (1|Com_nail) + (1|Com_muc) + (1|Com_scalp) + Com_BMI +
    (1|Hty_fh) + (1|Hty_exanth_1.0) + (1|Hty_exanth_2.0) +(1|Hty_self) + (1|Hty_photos) +
    (1|Hty_progr) + (1|Hty_smoker) + Hty_sport + Hty_alc + (1|Hty_onset) + Hty_pruritus +
    Hty_chronic)
varPart.res0.5 <- variancePartition::fitExtractVarPartModel(
  exprObj=v.object, formula=formula.design, data=info )


# 2.3.2 res 0.8
formula.design = formula(
  ~ (1|Molecular.subtype.res0.8) + (1|diag) + (1|batchID) + 
    (1|Pattern)  +  age.x + (1|Sex.x) +
    Morph_erythema + Morph_Elev + Morph_papu + (1|Morph_vesicles) + Morph_pustules +
    Morph_scaling + Morph_dryness + (1|Diag_PGA.x) + Diag_PASI + (1|Diag_SCORAD) + 
    Diag_DLQI + (1|Therapeutic.target_IL.17) + (1|Therapeutic.target_IL.23) + 
    (1|Therapeutic.target_MTX) + (1|Therapeutic.target_TNF) + The_dor + The_sor + 
    The_leuco + The_Hb + The_lympho + The_weeks + The_granulo + The_eosino + The_crea +
    The_GPT + (1|Lab_csa) + (1|Lab_specificIgE) + Lab_leuco + Lab_Hb + Lab_lympho +
    Lab_granulo + Lab_eosino + Lab_crea + Lab_GPT + Lab_totalIgE +  Hist_SE + 
    Hist_Para_quali + Hist_Aka_quali + Hist_Granu + Hist_Serum + Hist_Dist_Lympho +
    Hist_Microabscess + (1|Hist_Bacteria) + Hist_Neutro_quali + (1|Hist_Keratoses) +
    Hist_ID_distribution + Hist_ID_quali + Hist_Nr_Keratoses + Hist_Hyper_quant +
    Hist_Para_quant + Hist_Ortho_quant + Hist_Aka_quant + Hist_Spongiosis + 
    Hist_Cap + Hist_Lympho + Hist_Exocytosis + Hist_Neutro_quant + Hist_Eosino +
    (1|Hist_ID_quant) + (1|Hist_LCV) + Hist_Mucin + (1|Com_arths) + (1|Com_ah) + 
    (1|Com_iaa) + (1|Com_rca) + (1|Com_dm) + (1|Com_ast) + (1|Com_hepa) + (1|Com_renal) + 
    (1|Com_tonsil) + (1|Com_ibd) + (1|Com_nail) + (1|Com_muc) + (1|Com_scalp) + Com_BMI +
    (1|Hty_fh) + (1|Hty_exanth_1.0) + (1|Hty_exanth_2.0) +(1|Hty_self) + (1|Hty_photos) +
    (1|Hty_progr) + (1|Hty_smoker) + Hty_sport + Hty_alc + (1|Hty_onset) + Hty_pruritus +
    Hty_chronic)
varPart.res0.8 <- variancePartition::fitExtractVarPartModel(
  exprObj=v.object, formula=formula.design, data=info )


# 2.3.3 res 0.9
formula.design = formula(
  ~ (1|Molecular.subtype.res0.9) + (1|diag) + (1|batchID) + 
    (1|Pattern)  +  age.x + (1|Sex.x) +
    Morph_erythema + Morph_Elev + Morph_papu + (1|Morph_vesicles) + Morph_pustules +
    Morph_scaling + Morph_dryness + (1|Diag_PGA.x) + Diag_PASI + (1|Diag_SCORAD) + 
    Diag_DLQI + (1|Therapeutic.target_IL.17) + (1|Therapeutic.target_IL.23) + 
    (1|Therapeutic.target_MTX) + (1|Therapeutic.target_TNF) + The_dor + The_sor + 
    The_leuco + The_Hb + The_lympho + The_weeks + The_granulo + The_eosino + The_crea +
    The_GPT + (1|Lab_csa) + (1|Lab_specificIgE) + Lab_leuco + Lab_Hb + Lab_lympho +
    Lab_granulo + Lab_eosino + Lab_crea + Lab_GPT + Lab_totalIgE +  Hist_SE + 
    Hist_Para_quali + Hist_Aka_quali + Hist_Granu + Hist_Serum + Hist_Dist_Lympho +
    Hist_Microabscess + (1|Hist_Bacteria) + Hist_Neutro_quali + (1|Hist_Keratoses) +
    Hist_ID_distribution + Hist_ID_quali + Hist_Nr_Keratoses + Hist_Hyper_quant +
    Hist_Para_quant + Hist_Ortho_quant + Hist_Aka_quant + Hist_Spongiosis + 
    Hist_Cap + Hist_Lympho + Hist_Exocytosis + Hist_Neutro_quant + Hist_Eosino +
    (1|Hist_ID_quant) + (1|Hist_LCV) + Hist_Mucin + (1|Com_arths) + (1|Com_ah) + 
    (1|Com_iaa) + (1|Com_rca) + (1|Com_dm) + (1|Com_ast) + (1|Com_hepa) + (1|Com_renal) + 
    (1|Com_tonsil) + (1|Com_ibd) + (1|Com_nail) + (1|Com_muc) + (1|Com_scalp) + Com_BMI +
    (1|Hty_fh) + (1|Hty_exanth_1.0) + (1|Hty_exanth_2.0) +(1|Hty_self) + (1|Hty_photos) +
    (1|Hty_progr) + (1|Hty_smoker) + Hty_sport + Hty_alc + (1|Hty_onset) + Hty_pruritus +
    Hty_chronic)
varPart.res0.9 <- variancePartition::fitExtractVarPartModel(
  exprObj=v.object, formula=formula.design, data=info )

# TODO @Guelce: please save the varPart object so that I can use it for further analysis :)

# vp <- sortCols( varPart )
# plotPercentBars( vp[1:10,] )
# plotVarPart( vp )


######################################################################
#                     Removing batch effects                         #
######################################################################
# 3. Removing batch effects prior to variancePartition to assess the
# influence of biological variables better, subtract out effect of technical factors
# with linear mixed model
# Memory efficient: extract residuals directly without storing intermediate results
residList <- fitVarPartModel( exprObj=v.object, formula=~ (1|batchID) + 
                                age.x + (1|Sex.x), data=info )
# convert list to matrix
residMatrix = do.call(rbind, residList)

# fit model on residuals
form.v2 <- formula(
  ~ (1|Molecular.subtype.res0.5) + (1|diag) + (1|Pattern) +
    Morph_erythema + Morph_Elev + Morph_papu + (1|Morph_vesicles) + Morph_pustules +
    Morph_scaling + Morph_dryness + (1|Diag_PGA.x) + Diag_PASI + (1|Diag_SCORAD) + 
    Diag_DLQI + (1|Therapeutic.target_IL.17) + (1|Therapeutic.target_IL.23) + 
    (1|Therapeutic.target_MTX) + (1|Therapeutic.target_TNF) + The_dor + The_sor + 
    The_leuco + The_Hb + The_lympho + The_weeks + The_granulo + The_eosino + The_crea +
    The_GPT + (1|Lab_csa) + (1|Lab_specificIgE) + Lab_leuco + Lab_Hb + Lab_lympho +
    Lab_granulo + Lab_eosino + Lab_crea + Lab_GPT + Lab_totalIgE +  Hist_SE + 
    Hist_Para_quali + Hist_Aka_quali + Hist_Granu + Hist_Serum + Hist_Dist_Lympho +
    Hist_Microabscess + (1|Hist_Bacteria) + Hist_Neutro_quali + (1|Hist_Keratoses) +
    Hist_ID_distribution + Hist_ID_quali + Hist_Nr_Keratoses + Hist_Hyper_quant +
    Hist_Para_quant + Hist_Ortho_quant + Hist_Aka_quant + Hist_Spongiosis + 
    Hist_Cap + Hist_Lympho + Hist_Exocytosis + Hist_Neutro_quant + Hist_Eosino +
    (1|Hist_ID_quant) + (1|Hist_LCV) + Hist_Mucin + (1|Com_arths) + (1|Com_ah) + 
    (1|Com_iaa) + (1|Com_rca) + (1|Com_dm) + (1|Com_ast) + (1|Com_hepa) + (1|Com_renal) + 
    (1|Com_tonsil) + (1|Com_ibd) + (1|Com_nail) + (1|Com_muc) + (1|Com_scalp) + Com_BMI +
    (1|Hty_fh) + (1|Hty_exanth_1.0) + (1|Hty_exanth_2.0) +(1|Hty_self) + (1|Hty_photos) +
    (1|Hty_progr) + (1|Hty_smoker) + Hty_sport + Hty_alc + (1|Hty_onset) + Hty_pruritus +
    Hty_chronic)
varPartResid.res0.5 <- variancePartition::fitExtractVarPartModel(
  exprObj=residMatrix, formula=form.v2, data=info )

form.v2 <- formula(
  ~ (1|Molecular.subtype.res0.8) + (1|diag) + (1|Pattern) +
    Morph_erythema + Morph_Elev + Morph_papu + (1|Morph_vesicles) + Morph_pustules +
    Morph_scaling + Morph_dryness + (1|Diag_PGA.x) + Diag_PASI + (1|Diag_SCORAD) + 
    Diag_DLQI + (1|Therapeutic.target_IL.17) + (1|Therapeutic.target_IL.23) + 
    (1|Therapeutic.target_MTX) + (1|Therapeutic.target_TNF) + The_dor + The_sor + 
    The_leuco + The_Hb + The_lympho + The_weeks + The_granulo + The_eosino + The_crea +
    The_GPT + (1|Lab_csa) + (1|Lab_specificIgE) + Lab_leuco + Lab_Hb + Lab_lympho +
    Lab_granulo + Lab_eosino + Lab_crea + Lab_GPT + Lab_totalIgE +  Hist_SE + 
    Hist_Para_quali + Hist_Aka_quali + Hist_Granu + Hist_Serum + Hist_Dist_Lympho +
    Hist_Microabscess + (1|Hist_Bacteria) + Hist_Neutro_quali + (1|Hist_Keratoses) +
    Hist_ID_distribution + Hist_ID_quali + Hist_Nr_Keratoses + Hist_Hyper_quant +
    Hist_Para_quant + Hist_Ortho_quant + Hist_Aka_quant + Hist_Spongiosis + 
    Hist_Cap + Hist_Lympho + Hist_Exocytosis + Hist_Neutro_quant + Hist_Eosino +
    (1|Hist_ID_quant) + (1|Hist_LCV) + Hist_Mucin + (1|Com_arths) + (1|Com_ah) + 
    (1|Com_iaa) + (1|Com_rca) + (1|Com_dm) + (1|Com_ast) + (1|Com_hepa) + (1|Com_renal) + 
    (1|Com_tonsil) + (1|Com_ibd) + (1|Com_nail) + (1|Com_muc) + (1|Com_scalp) + Com_BMI +
    (1|Hty_fh) + (1|Hty_exanth_1.0) + (1|Hty_exanth_2.0) +(1|Hty_self) + (1|Hty_photos) +
    (1|Hty_progr) + (1|Hty_smoker) + Hty_sport + Hty_alc + (1|Hty_onset) + Hty_pruritus +
    Hty_chronic)
varPartResid.res0.8 <- variancePartition::fitExtractVarPartModel(
  exprObj=residMatrix, formula=form.v2, data=info )

form.v2 <- formula(
  ~ (1|Molecular.subtype.res0.9) + (1|diag) + (1|Pattern) +
    Morph_erythema + Morph_Elev + Morph_papu + (1|Morph_vesicles) + Morph_pustules +
    Morph_scaling + Morph_dryness + (1|Diag_PGA.x) + Diag_PASI + (1|Diag_SCORAD) + 
    Diag_DLQI + (1|Therapeutic.target_IL.17) + (1|Therapeutic.target_IL.23) + 
    (1|Therapeutic.target_MTX) + (1|Therapeutic.target_TNF) + The_dor + The_sor + 
    The_leuco + The_Hb + The_lympho + The_weeks + The_granulo + The_eosino + The_crea +
    The_GPT + (1|Lab_csa) + (1|Lab_specificIgE) + Lab_leuco + Lab_Hb + Lab_lympho +
    Lab_granulo + Lab_eosino + Lab_crea + Lab_GPT + Lab_totalIgE +  Hist_SE + 
    Hist_Para_quali + Hist_Aka_quali + Hist_Granu + Hist_Serum + Hist_Dist_Lympho +
    Hist_Microabscess + (1|Hist_Bacteria) + Hist_Neutro_quali + (1|Hist_Keratoses) +
    Hist_ID_distribution + Hist_ID_quali + Hist_Nr_Keratoses + Hist_Hyper_quant +
    Hist_Para_quant + Hist_Ortho_quant + Hist_Aka_quant + Hist_Spongiosis + 
    Hist_Cap + Hist_Lympho + Hist_Exocytosis + Hist_Neutro_quant + Hist_Eosino +
    (1|Hist_ID_quant) + (1|Hist_LCV) + Hist_Mucin + (1|Com_arths) + (1|Com_ah) + 
    (1|Com_iaa) + (1|Com_rca) + (1|Com_dm) + (1|Com_ast) + (1|Com_hepa) + (1|Com_renal) + 
    (1|Com_tonsil) + (1|Com_ibd) + (1|Com_nail) + (1|Com_muc) + (1|Com_scalp) + Com_BMI +
    (1|Hty_fh) + (1|Hty_exanth_1.0) + (1|Hty_exanth_2.0) +(1|Hty_self) + (1|Hty_photos) +
    (1|Hty_progr) + (1|Hty_smoker) + Hty_sport + Hty_alc + (1|Hty_onset) + Hty_pruritus +
    Hty_chronic)
varPartResid.res0.9 <- variancePartition::fitExtractVarPartModel(
  exprObj=residMatrix, formula=form.v2, data=info )

# TODO @Guelce: please save also the varPartResid object so that I can use it for further analysis :)


