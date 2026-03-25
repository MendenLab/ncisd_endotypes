rm(list = ls(all=TRUE))

library(edgeR)
library(limma)
library(DESeq2)
library('xlsx')
library(dplyr)  
library(ggsankey)
library('ggalluvial')
library(tidyr)
library('alluvial')  # https://github.com/mbojan/alluvial
library('forcats')


# To estimate the variance in lesional gene expression explained by differences in 
# diagnosis, gender, age and other clinical attributes, 
# the RNAseq count data was first transformed using voom transformation from R package Limma and then modelled with a linear mixed effect model using bioconductor package variancePartition.


#' @description Load .rds object containing lesion and non-lesion samples including metaData
#' @param input.dir Path to count matrix and meta.data
def.load_data <- function(input.dir, filename='dds_highQual_Sexandbatchcorrected_v04.rds') 
{
  
  # Load raw expression data but filtered for protein transcribing genes and TPM > 1
  # Design: ~ SampleType + Sex.x + batch
  dds <- readRDS(file.path(input.dir, 'count_matrix', 'data_freeze', 'preprocessed', 'L_and_NL',
                           filename))
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
output.dir <- file.path(output.dir, 'SankeyDiagramm', 'Figure_2C', Sys.Date())
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

# 1. Data preparation
obs.name <- "healthysamp_diag_v2"
# 1.1 Load data
dds <- def.load_data(input.dir="/Volumes/CH__data/Projects/Eyerich_AG_projects/BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer", 
                     filename='dds_highQual_Sexandbatchcorrected__Endotypes_renamed_v01.rds')

# Load therapy response data
df.metadata <- read.csv('/Users/christina.hillig/PycharmProjects/Eyerich_ERCGrant/Molecular_subtypes/output/Figure_UMAP_Hoverover/2023-03-03/MetaData.csv')
df.metadata <- df.metadata[, colSums(is.na(df.metadata)) != nrow(df.metadata)]
ind.muc.dds <- match(df.metadata$Helmholtz_identifyer, colnames(dds[ , dds$sampleType == 'D']))
# Add therapy information to dds object
for (col.therapy in colnames(df.metadata)[seq(110, 137, 2)]) 
{
  dds[[col.therapy]] <- ""
  dds[[col.therapy]][dds$sampleType == 'D'][ind.muc.dds] <- df.metadata[[col.therapy]]
}


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


###############################################################
#               Diag - Pattern - Molecular_subtype            #
###############################################################
# Create dataframe for sankey diagramm 
df.sankey <- y.object$samples %>% dplyr::select(Molecular_subtype, sampleType, diag, Pattern)
df.sankey$MUCIDs = row.names(y.object$samples)
# Select only lesion samples
df.sankey <- df.sankey[df.sankey$sampleType == 'D', ]
df.sankey$sampleType <- NULL

df.sankey$Molecular_subtype <- factor(df.sankey$Molecular_subtype, 
                         levels = c(
                           'E1', 'E2', 'E3', 'E4', 'E5', 'E6', 'E7', 'E8', 
                           'E9', 'E10', 'E11', 'E12', 'E13'))
# Get names
Molecular_subtype.names <- levels(as.factor(df.sankey$Molecular_subtype))

# Refactor pattern
df.sankey$Pattern <- as.factor(as.character(df.sankey$Pattern))
pattern.names <- levels(df.sankey$Pattern)

# Reorder diag and use abbreviations
df.sankey$diag <- as.factor(as.character(df.sankey$diag))
df.sankey$diag <- factor(df.sankey$diag, 
                         levels = c(
                           'lichen planus', 'lupus erythematosus', 'lichenoid drug reaction', 
                           'eczema', 'prurigo simplex subacuta', 'bullous pemphigoid', 
                           'psoriasis','pityriasis rubra pilaris', 'morphea', 
                           'venous ulcer', 'systemic sclerosis', 'granuloma annulare',
                           'sarcoidosis', 'psoriasis pustulosa', 'pyoderma gangrenosum', 
                           'cutaneous lymphoma', 'cutaneous side effects of biologics', 'darier disease', 
                           'keratosis lichenoides chronica', 'erythrodermia', 'parapsoriasis', 'undefined'))
nm1 <- setNames(c('LP', 'LE', 'LDR', 'Eczema', 'PSS', 'BP', 'Pso', 'PRP', 'Morphea', 'VU', 'SS', 'GA', 
                  'Sarco.', 'PP', 'PG', 'CTL', 'CSEoB', 'DD', 'KLC', 'Erythro.', 'ParaPso', 'UD'),
                c('lichen planus', 'lupus erythematosus', 'lichenoid drug reaction', 'eczema', 'prurigo simplex subacuta', 
                  'bullous pemphigoid', 'psoriasis','pityriasis rubra pilaris', 'morphea', 
                  'venous ulcer', 'systemic sclerosis', 'granuloma annulare',
                  'sarcoidosis', 'psoriasis pustulosa', 'pyoderma gangrenosum', 'cutaneous lymphoma',
                  'cutaneous side effects of biologics', 'darier disease', 
                  'keratosis lichenoides chronica', 'erythrodermia', 'parapsoriasis', 'undefined'))
df.sankey$diag <- unname(nm1[df.sankey$diag])
df.sankey$diag <- factor(
  df.sankey$diag,  levels = c('LP', 'LE', 'LDR', 'Eczema', 'PSS', 'BP', 'Pso', 'PRP', 'Morphea', 'VU', 'SS', 'GA', 
                              'Sarco.', 'PP', 'PG', 'CTL', 'CSEoB', 'DD', 'KLC', 'Erythro.', 'ParaPso', 
                              'UD'))
diag.names <- levels(df.sankey$diag)

df.sankey.mutated <- df.sankey %>% 
  count(diag, Pattern, Molecular_subtype) %>% 
  mutate(diag = fct_relevel(as.factor(diag), rev(diag.names)))  %>%
  mutate(Pattern = fct_relevel(as.factor(Pattern), rev(pattern.names)))  %>%
  mutate(Molecular_subtype = fct_relevel(as.factor(Molecular_subtype), 
                                                rev(Molecular_subtype.names)))

pdf(file = file.path(output.dir, 
                     paste0('SankeyDiagramm_Diag_Pattern_Molecular_subtype', '.pdf')), 
    width = 10, height = 14)
print(alluvial::alluvial(
  df.sankey.mutated %>% select(-n),
  freq = df.sankey.mutated$n, border = NA, alpha = 0.5,
  col=case_when(df.sankey.mutated$Molecular_subtype == "E1" ~ "#006ddb",
                df.sankey.mutated$Molecular_subtype == "E2" ~ "#b6dbff",
                df.sankey.mutated$Molecular_subtype == "E3" ~ "#004949",
                df.sankey.mutated$Molecular_subtype == "E4" ~ "#009292",
                df.sankey.mutated$Molecular_subtype == "E5" ~ "#ff6db6",
                df.sankey.mutated$Molecular_subtype == "E6" ~ "#490092",
                df.sankey.mutated$Molecular_subtype == "E7" ~ "#b66dff",
                df.sankey.mutated$Molecular_subtype == "E8" ~ "#000000",
                df.sankey.mutated$Molecular_subtype == "E9" ~ "#920000",
                df.sankey.mutated$Molecular_subtype == "E10" ~ "#E69F00",
                df.sankey.mutated$Molecular_subtype == "E11" ~ "#D55E00",
                df.sankey.mutated$Molecular_subtype == "E12" ~ "#8B4513",
                df.sankey.mutated$Molecular_subtype == "E13" ~ "#999999",
                TRUE ~ "black"),
  cex = 0.75,
  axis_labels = c("Diagnosis", "Pattern", "Endotypes"))) 
dev.off()




