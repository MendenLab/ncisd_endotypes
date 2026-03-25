library(tibble)

do_add.column <- function(df, entrezid, gene.name, ensembl.id) 
{
  # add genenames and entrezid to dataframe
  df <- tibble::add_column(df, gene_name = ensembl.id, .after = "hgnc_symbol")
  df <- tibble::add_column(df, entrezid = entrezid, .after = "hgnc_symbol")
  df <- tibble::add_column(df, gene_name = gene.name, .after = "hgnc_symbol")
  
  return(df)
}


do.batch_correction <- function(func_dds, batch.effect, biological.corvar, 
                                color.obs, shape.obs, color.col, shape.col, func_save.dir) 
{
  # https://rnabio.org/module-03-expression/0003/05/01/Batch-Correction/ 
  uncorrected_data <- counts(func_dds)
  # drop levels 
  func_dds$Donor <- droplevels(func_dds$Donor)
  func_dds$Stimuli <- droplevels(func_dds$Stimuli)
  
  # TODO apply batch correction 
  # include condition (group variable)
  # Raw count matrix from genomic studies (dimensions gene x sample)
  # BiocManager::install("sva")
  corrected_data <- sva::ComBat_seq(
    uncorrected_data, batch=batch.effect, 
    covar_mod=colData(func_dds)[, biological.corvar], group=NULL, full_mod=TRUE)
  
  
  #compare dimensions of corrected and uncorrected data sets
  print(dim(uncorrected_data))
  print(dim(corrected_data))
  
  #visually compare values of corrected and uncorrected data sets
  print(head(uncorrected_data))
  print(head(corrected_data))
  
  # create new DESeq2 object with corrected counts
  # plotPCA by default uses the top 500 most variable genes prior to prcomp
  pca_corrected_obj.allgenes <- stats::prcomp(t(corrected_data))
  print(head(pca_corrected_obj.allgenes))
  # pull PCA values out of the PCA object
  percentVar <- pca_corrected_obj.allgenes$sdev^2 / sum(pca_corrected_obj.allgenes$sdev^2)
  intgroup.df <- as.data.frame(colData(func_dds)[ , c(color.col, shape.col)])
  # group <- if (length(c(color.col, shape.col)) > 1) {
  #   factor(apply(intgroup.df, 1, paste, collapse = ":"))}
  pca_corrected <- data.frame(PC1 = pca_corrected_obj.allgenes$x[, 1], 
                              PC2 = pca_corrected_obj.allgenes$x[, 2], 
                              name = colnames(func_dds))
  pca_corrected <- cbind(pca_corrected, intgroup.df)
  
  return(list(corrected_data, pca_corrected))
}


def.prepare_subtypes <- function(df.file, obs='Molecular.subtypes') 
{
  # Save subtypes for seed 0
  subtypes.list <- df.file[1, obs]
  subtypes.list <- (gsub("\\[|\\]|[\r\n]|'|,", "", subtypes.list))
  subtypes.list <-  strsplit(subtypes.list, " ")[[1]]
  
  return(subtypes.list)
}


add.endotypes_LNL <- function(dds, df) 
{
  list.subtypes <- def.prepare_subtypes(df.file = df)
  # To start at 1 instead of 0 - might be removed in future
  list.subtypes = as.integer(list.subtypes) + 1
  # drop NAs
  list.subtypes <- list.subtypes[!is.na(list.subtypes)]
  
  # Remove samples which are not part of the lesion dds object ..
  mucids.L <- def.prepare_subtypes(df.file=df, obs='MUC.IDs') 
  mucids.NL <- dds[,(dds$sampleType == 'H')]$sampleID
  mucids <- c(mucids.L, mucids.NL)
  dds <- dds[ , dds$sampleID %in% mucids]
  # # Find matching MUC IDs
  ind.dds <- match(mucids.L, colnames(dds[ , dds$sampleType == 'D']))
  
  # Add subtypes
  dds$Molecular_subtype <- 'non-lesional' 
  dds$Molecular_subtype[dds$sampleType == 'D'][ind.dds] <- list.subtypes
  dds$Molecular_subtype <- as.factor(dds$Molecular_subtype)
  
  
  return(dds)
  
}

add.endotypes_L <- function(dds, df) 
{
  list.subtypes <- def.prepare_subtypes(df.file = df)
  # To start at 1 instead of 0 - might be removed in future
  list.subtypes = as.integer(list.subtypes) + 1
  # drop NAs
  list.subtypes <- list.subtypes[!is.na(list.subtypes)]
  
  # Remove samples which are not in the only lesion object ..
  mucids.L <- def.prepare_subtypes(df.file=df, obs='MUC.IDs') 
  dds <- dds[ , dds$sampleID %in% mucids.L]
  # # Find matching MUC IDs
  ind.dds <- match(mucids.L, colnames(dds[ , dds$sampleType == 'D']))
  
  # Add Endotypes
  dds$`Molecular_subtype` <- 'L'
  dds$`Molecular_subtype`[dds$sampleType == 'D'][ind.dds] <- list.subtypes
  dds$`Molecular_subtype` <- as.factor(dds$`Molecular_subtype`)
  
  return(dds)
}


drop_lesion_samples_wo_NL_samples <- function(dds, obs.name) 
{
  names.keep.vector <- c()
  for (subtype in levels(dds[[obs.name]]))
  {
    mask.subtype <- dds[[obs.name]] == subtype
    names.keep <- names(which(table(dds$PatientID[mask.subtype]) >= 2))
    
    names.keep.vector <- c(names.keep.vector, names.keep)
  }
  dds <- dds[, dds$PatientID %in% names.keep.vector]
  dds$PatientID <- base::droplevels(dds$PatientID)
  
  return(dds)
}


def.rename_nonlesion_samples_to_paired_lesion_name <- function(dds, ref.obs.name, obs.name, group.name) 
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
    # Find healthy partner via Patient ID == Pseudo.ID
    patient.id <- dds$PatientID[mask.observable]
    # multiple matching using which and %in%, as match only return first encounter 
    mask.patients <- which(dds$PatientID %in% patient.id)
    # TODO remove those which are the patient but are clustered in a different cluster  
    mask.patients <- mask.patients[
      (dds[[group.name]][mask.patients] == observable) | (
        dds[[group.name]][mask.patients] == 'non-lesional')]
    dds[[obs.name]][mask.patients] <- observable
  }
  dds[[obs.name]] <- as.factor(dds[[obs.name]])
  
  return(dds)
}
