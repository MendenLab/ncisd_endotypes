# Function of script: Prepare 2D/3D count matrix
# R version 4.1.1 (2021-08-10)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Catalina 10.15.7

# Load libraries
library(biomaRt)
library(org.Hs.eg.db) 
library(dplyr)
library(tidyr)

general.dir <- '/Volumes/CH__data/Projects/Eyerich_AG_projects/BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer'
source(file.path(general.dir, 'analysis/Molecular_subtypes/r_scripts/DGE_analysis/utils.R'))


get_biotypes <- function(names.genes, filter.option) 
{
  ## set up connection to ensembl database and specify a data set to use
  ensembl <- biomaRt::useMart(
    biomart = "ENSEMBL_MART_ENSEMBL",  # or "genes"??
    dataset="hsapiens_gene_ensembl",
    host = 'https://may2021.archive.ensembl.org' # version 104
  )
  
  #find gene attributes
  bm <- getBM(attributes = 
                c("hgnc_symbol", # HGNC approves gene names (HUGO-gene names)
                  "ensembl_gene_id",
                  'entrezgene_id',
                  "external_gene_name",
                  "description",
                  "gene_biotype",
                  "chromosome_name"),
              filter= filter.option,# "external_gene_name", # "ensembl_gene_id",
              values= names.genes,
              mart=ensembl, uniqueRows=TRUE)
  
  # Source: https://www.biostars.org/p/447677/
  # The fundamental reason for the mapping of one HGNC gene symbol to many Ensembl genes is 
  # the mismatch between gene definitions in HGNC and Ensembl. Gene definition in Ensembl is
  # locus-based because it is associated with annotation of a reference genome whereas HGNC 
  # has this definition:
  # --> A gene is defined as: "a DNA segment that contributes to phenotype/function. In the
  # absence of demonstrated function a gene may be characterized by sequence, transcription 
  # or homology".
  # ==> Note that the HGNC definition can apply to multiple loci in the genome hence to 
  # multiple Ensembl genes.
  
  
  # arrange them for their gene ids
  bm <- arrange(bm, hgnc_symbol)
  
  # compare length of gene names in received data against biomart results
  print(length(names.genes))
  print(length(bm$hgnc_symbol))
  
  return(bm)
}


map_genesymbol_entrezid_genename <- function(gene.symbols) 
{
  #'
  #' @description 
  #' @param 
  #' @return 
  
  genename <- vector()
  entrezid <- vector()
  ensembl.id <- vector()
  for (sym_c in 1:length(gene.symbols)) 
  {
    # print('***')
    # print(gene.symbols[sym_c])
    entrezid[sym_c] <- as.character(mget(gene.symbols[sym_c], org.Hs.eg.db::org.Hs.egSYMBOL2EG,
                                         ifnotfound = NA))
    # print(entrezid[sym_c])
    genename[sym_c] <- as.character(mget(entrezid[sym_c], org.Hs.eg.db::org.Hs.egGENENAME,
                                         ifnotfound = NA))
    
    # convert entrez id to ensembl
    ensembl.id[sym_c] <- as.character(mget(entrezid[sym_c],
                                           org.Hs.eg.db::org.Hs.egENSEMBL,
                                           ifnotfound = NA)) 
  }
  
  
  return(list(entrezid, genename, ensembl.id))
}


add__biotypes <- function(dge.object, bio.types) 
{
  func_raw_counts <- counts(dge.object)
  
  # make sure rows of bio.types line up with the countmatrix of the dge.object
  if (all(length(rownames(dge.object)) == length(bio.types$hgnc_symbol)))
  {
    #create a new datarframe with a column of the gene ids, in the order of the countmatrix
    annot <- data.frame(hgnc_symbol = rownames(func_raw_counts), 
                        row.names = rownames(func_raw_counts))
    annot$hgnc_symbol <- as.character(annot$hgnc_symbol)
    
    # join dataframes 
    # by ensures that each hgnc_symbol in annot gets its corresponding attributes in bio.types
    gene_attributes <- dplyr::left_join(annot, bio.types, by = "hgnc_symbol")
    ### add the annotation to the DESeq2 table dge.object
    mcols(rowRanges(dge.object)) <- gene_attributes
    
    
    #extract biotype + counts in one element
    counts_and_attributes <- cbind(mcols(rowRanges(dge.object))$gene_biotype, assay(dge.object))
    counts_and_attributes <- as.data.frame(counts_and_attributes)
    
    ggplot(data = counts_and_attributes, aes(x=biotype, y=allcounts, color=biotype)) +
      geom_boxplot() +
      scale_y_continuous("counts") +
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank())
    
  } else 
  {
    # else convert gene names to entriez id and esembl id via library(org.Hs.eg.db) 
    entrez.name.ensembl <- map_genesymbol_entrezid_genename(gene.symbols = rownames(dge.object)) 
    
    #create a new dataframe with a column of the gene ids, in the order of the countmatrix
    annot <- data.frame(hgnc_symbol = rownames(func_raw_counts))
    annot$hgnc_symbol <- as.character(annot$hgnc_symbol)
    
    annot <- do_add.column(df = annot, entrezid = entrez.name.ensembl[[1]],
                           gene.name = entrez.name.ensembl[[2]], 
                           ensembl.id = entrez.name.ensembl[[3]]) 
    
    # rename columns
    annot <- annot %>% rename('description' = 'gene_name.1', 'ENSEMBL' = 'gene_name')
    # annot length: 24538
    
    # get biotypes
    func_biotypes.list = get_biotypes(names.genes = entrez.name.ensembl[[3]],
                                      filter.option = 'ensembl_gene_id')
    
    # keep only protein coding genes
    func_biotypes.list <- func_biotypes.list[
      func_biotypes.list$gene_biotype == 'protein_coding', ]
    index_biolist <- match(func_biotypes.list$hgnc_symbol, annot$hgnc_symbol)
    index_biolist <- index_biolist[!is.na(index_biolist)]
    annot <- annot[index_biolist, ]
    func_biotypes.list <- func_biotypes.list[index_biolist, ]
    
    # combine dataframes
    # by ensures that each hgnc_symbol in annot gets its 
    # corresponding attributes in func_biotypes.list
    annot <- tibble::add_column(
      annot, gene_name = func_biotypes.list$gene_biotype, .after = "hgnc_symbol")
    
    ### add the annotation to the DESeq2 table dds
    # all(index_biolist == index_dge.object) -> TRUE
    index_dge.object <- match(annot$hgnc_symbol, row.names(dge.object))
    # subset dge.object
    dge.object <- dge.object[index_dge.object]
    rowData(dge.object) <- annot
  }
  
  return(dge.object)
}



def.prepare_metadata <- function(func_meta.data, func_count.matrix) 
{
  func_filtered.count.matrix <- func_count.matrix[ , func_meta.data$Helmholtz.identifyer]
  
  print("Same order of Samples names in metaData and Count matrix: ")
  print(all(func_meta.data$Helmholtz.identifyer == colnames(func_filtered.count.matrix)))
  
  return(list(func_meta.data, func_filtered.count.matrix))
}


filter_zerocounts <- function(data)
{
  # input: DESeq2 object
  # minimal filtering rule: removing rows that have no counts, or only a single count across all samples
  print(nrow(data))
  # only keep genes with added counts >= 1
  keep <- rowSums(counts(data)) >= 1
  data <- data[keep, ]
  print(nrow(data))
  return(data)
}


TPM_filter <-  function(Raw_Count_Matrix, Ideal_Genelength = NULL, minsumTPM_pergene) {
  # Input is dds object -> read out counts
  raw_counts <- as.data.frame(counts(Raw_Count_Matrix))
  row.names(raw_counts) <- rowData(Raw_Count_Matrix)$ENSEMBL
  
  if (missing(minsumTPM_pergene)) {
    minsumTPM_pergene<- ncol(raw_counts)
    # tpm<1 (transcripts <1) in "minTPM_pergene" samples
  }
  raw_counts$ensembl_id <- rownames(raw_counts)
  # remove .xx from ENSEMBL ID
  raw_counts$ensembl_id <- gsub("\\.", "_", raw_counts$ensembl_id)
  raw_counts$ensembl_id <- sub("_.*", "", raw_counts$ensembl_id)
  # remove duplicated rows
  raw_counts <- raw_counts[!duplicated(raw_counts$ensembl_id), ]
  
  if (is.null(Ideal_Genelength)) {
    Ideal_Genelength <- readRDS("/Users/christina.hillig/R_studio/Ideal_Genelength.rds")
    Final_Count_Matrix <- dplyr::inner_join(raw_counts, Ideal_Genelength, by= "ensembl_id")
  } else { # a idealgenelength dataframe is provided
    rownames(Ideal_Genelength) <- c("ensembl_id", "ideal_length")
    Final_Count_Matrix <- dplyr::inner_join(raw_counts, Ideal_Genelength, by= "ensembl_id")
  }
  
  tpm_tibble <- Final_Count_Matrix %>% tidyr::pivot_longer(
    cols = starts_with("M"), names_to = "sampleID", values_to = "count")  %>%
    group_by(sampleID) %>% mutate(
      rpk = count / (ideal_length / 1000), libSize = sum(rpk), tpm = rpk / libSize *1e6 ) %>%
    mutate(minTPM= min(tpm[tpm > 0]), log10TPM = log10(tpm + minTPM))
  
  tpm_tibble2 <- tpm_tibble %>% ungroup() %>% group_by(ensembl_id) %>%
    dplyr::filter( sum(tpm < 1) < minsumTPM_pergene ) ## filter genes which have tpm < 1 , i.e. which have even 1 transcript absent in all samples
  
  cnt.tib <- tpm_tibble2 %>% dplyr::select(ensembl_id, sampleID, count) %>%
    tidyr::pivot_wider(names_from = sampleID, values_from = count )
  
  Final_Count_Matrix <- cnt.tib %>% as.data.frame() %>% dplyr::select(-ensembl_id)
  rownames(Final_Count_Matrix) <- cnt.tib$ensembl_id
  
  
  # sort metaData and count matrix like gene features
  ### add the annotation to the DESeq2 table dds
  index_dge.object <- match(row.names(Final_Count_Matrix), rowData(Raw_Count_Matrix)$ENSEMBL)
  # remove Nas
  index_dge.object <- index_dge.object[!is.na(index_dge.object)]
  # subset dge.object
  Raw_Count_Matrix <- Raw_Count_Matrix[index_dge.object]
  
  index_dge.object <- match(rowData(Raw_Count_Matrix)$ENSEMBL, row.names(Final_Count_Matrix))
  index_dge.object <- index_dge.object[!is.na(index_dge.object)]
  # -> TRUE
  print(all(row.names(Final_Count_Matrix)[index_dge.object] == 
              rowData(Raw_Count_Matrix)$ENSEMBL))
  
  print(paste0("Dimension after filtering out TPM < 1: ", nrow(Raw_Count_Matrix)))
  
  return(Raw_Count_Matrix)
}


def.assess_quality <- function(dds.object, title, func_save.dir) 
{
  # How many reads were sequenced for each sample ( = library sizes)?
  print(colSums(counts(dds.object)))
  
  par(mfrow=c(1, 1))
  func.p1 <- colSums(counts(dds.object)) %>% barplot
  def.plot.save(
    plot.to.plot = func.p1, 
    file.title = paste0('Quality_assessment_', title, '.png'), 
    save.dir = func_save.dir, width=6, height = 6)
}


#' @note There is the assumption that some genes are not changing across conditions! 
#' • Size factors should be around 1.
#' • Normalized counts are calculated via countsgeneX, sampleA/sizefactorsampleA
def.assess_sizefactors <- function(dds.object, title, func_save.dir) 
{
  par(mfrow=c(1, 1))
  func.p1 <- plot( sizeFactors(dds.object), colSums(counts(dds.object)), # assess them
                   ylab = "library sizes", xlab = "size factors", cex = .6 )
  def.plot.save(
    plot.to.plot = func.p1, 
    file.title = paste0('SizeFactors_assesment_', title, '.png'), 
    save.dir = func_save.dir, width=6, height = 6)
}


#' @note The read counts normalized for sequencing depth can be accessed via counts(..., normalized = TRUE). Check whether the normalization helped adjust global differences between the samples.
def.assess_normalisation <- function(dds.object, file.title, func_save.dir) 
{
  ## extracting normalized counts
  counts.sf_normalized <- counts(dds.object, normalized=TRUE)
  ## setting up the plotting layout
  png(filename=file.path(
    func_save.dir, paste0('SizeFactor_normed_counts_per_sample_', file.title, '.png')), 
    units="in", width=6, height=6, res=300)
  par(mfrow=c(1, 2))
  ## adding the boxplots
  boxplot(counts.sf_normalized, main = "SizeFactor normalized", cex = .6) 
  boxplot(counts(dds.object), main = "Read counts only", cex = .6)
  dev.off()
  
  # to plot the two box plots next to each other
  par(mfrow=c(1, 2)) 
  ## boxplots of non-normalized counts
  boxplot(log2(counts(dds.object) + 1), notch=TRUE, main = "Non-normalized read counts", 
          ylab="log2(read counts)", cex = .6)
  ## boxplots of size-factor normalized values
  func.p1 = boxplot(log2(counts(dds.object, normalize= TRUE) +1), notch=TRUE, 
                    main = "Size-factor-normalized read counts", ylab="log2(read counts)", cex = .6)
  def.plot.save(
    plot.to.plot = func.p1, 
    file.title = paste0('Log2_SizeFactor_normed_counts_per_sample_', file.title, '.png'), 
    save.dir = func_save.dir, width=6, height = 6)
}


#' @note Do this before and after VST 
#' @param dds.object : deseq2.object 
#' 
def.check_homoskedasticity <- function(dds.object, vst.object, title, func_save.dir)
{
  ## generate the base meanSdPlot using sequencing depth normalized log2(read counts)
  log.norm.counts <- log2(counts(dds.object, normalized=TRUE) + 1) 
  ## generate the plot
  msd_plot <- vsn::meanSdPlot(log.norm.counts,
                              ranks=FALSE, # show the data on the original scale
                              plot = FALSE)
  ## since vsn::meanSdPlot generates a ggplot2 object, this can be
  ## manipulated in the usual ways
  msd_plot$gg +
    ggtitle("Sequencing depth normalized log2(read counts)") + ylab("standard deviation")
  
  def.plot.save(
    plot.to.plot = msd_plot$gg , 
    file.title = paste0('mean_sd_Plot_normed_counts_', title, '.png'), 
    save.dir = func_save.dir, width=6, height = 6)
  
  # Mean vs sd of VST counts
  msd_plot <-  vsn::meanSdPlot(assay(vst.object), ranks=FALSE, plot = FALSE) 
  msd_plot$gg + ggtitle("VST transformation") 
  
  def.plot.save(
    plot.to.plot = msd_plot$gg , 
    file.title = paste0('mean_sd_Plot_vst_counts_', title, '.png'), 
    save.dir = func_save.dir, width=6, height = 6)
  
  
}


# Literature: 
# https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/08_practical.pdf
