rm(list = ls(all=TRUE))
set.seed(1)

# Add to library path
.libPaths("r_libs")

library("MuSiC")
library('rhdf5')
# library('zellkonverter')
library('SingleCellExperiment')
library(openxlsx)
# library(SummarizedExperiment)
library(scran)
library(scuttle) # alternative to scran
library(scater)

library(tidyverse)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(forcats)


convert_obs_contents_to_df <- function(obs_list) {
  result <- list()
  
  for (name in names(obs_list)) {
    entry <- obs_list[[name]]
    
    if (is.list(entry) && all(c("categories", "codes") %in% names(entry))) {
      # Categorical variable: convert codes to factor with category labels
      result[[name]] <- factor(
        entry$codes, levels = seq_along(entry$categories) - 1, 
        labels = entry$categories)
      
    } else if (is.data.frame(entry)) {
      # Data frame: add columns with prefix to avoid name clashes
      for (col in names(entry)) {
        result[[paste0(name, "_", col)]] <- entry[[col]]
      }
      
    } else if (is.vector(entry) || is.matrix(entry)) {
      # Other vector or matrix types (optional, adjust as needed)
      result[[name]] <- entry
      
    } else {
      # Fallback: try to coerce to vector
      result[[name]] <- as.vector(entry)
    }
  }
  
  # Combine everything into one data frame
  return(as.data.frame(result))
}



def_get_obs <- function(path.to.h5file)
{
  # List all items in /obs/
  obs_contents <- rhdf5::h5read(path.to.h5file, "/obs")
  n_rows <- length(h5read(path.to.h5file, "/obs/index"))  # number of cells
  
  obs_df <- convert_obs_contents_to_df(obs_contents)
  
  return(obs_df)
}


def_load_data <- function(projectpath, projectname, filename, metadata_filename) 
{
  path.to.h5file <- file.path(projectpath, projectname, 'data', filename)
  path.to.metadatafile <- file.path(
    projectpath, projectname, 'data', metadata_filename)
  
  # Read out observations
  # adata.obs <- read.csv(file=path.to.metadatafile)
  adata.obs <- def_get_obs(path.to.h5file)
  # adata.obs <- as.data.frame(rhdf5::h5read(path.to.h5file, "/obs/"))
  # make barcode names unique
  # unique.barcodes <- make.names(adata.obs[[1]], unique = TRUE)
  
  # Read out gene infos
  adata.vars <- rhdf5::h5read(path.to.h5file, "/var/")
  
  # Read out filtered, raw counts
  counts_data <- rhdf5::h5read(path.to.h5file, "/layers/counts/data")
  counts_indices <- rhdf5::h5read(path.to.h5file, "/layers/counts/indices")
  counts_indptr <- rhdf5::h5read(path.to.h5file, "/layers/counts/indptr")
  # Get matrix shape
  n_rows <- length(h5read(path.to.h5file, "/obs/index"))  # number of cells
  n_cols <- length(h5read(path.to.h5file, "/var/index"))  # number of genes
  
  # Reconstruct sparse matrix
  counts_sparse <- new("dgRMatrix",
                       p = as.integer(counts_indptr),
                       j = as.integer(counts_indices),
                       x = as.numeric(counts_data),
                       Dim = as.integer(c(n_rows, n_cols)))
  # add row and column names
  rownames(counts_sparse) <- h5read(path.to.h5file, "/obs/index")
  colnames(counts_sparse) <- h5read(path.to.h5file, "/var/index")
  
  rownames(adata.obs) <- rownames(counts_sparse) 
  
  # Otional: Convert to dense matrix if needed
  # counts_dense <- as.matrix(counts_sparse)
  # counts_sparse = NULL 
  
  counts_sparse <- Matrix::t(counts_sparse)

  nrow(counts_sparse) == nrow(as.data.frame(
    adata.vars, row.names = rownames(counts_sparse)))
  ncol(counts_sparse) == nrow(adata.obs)
  
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts_sparse),
    colData = adata.obs, rowData = as.data.frame(
      adata.vars, row.names = rownames(counts_sparse)))
  
  return(sce)
}


def_load_bulk_data <- function(projectpath, projectname, filename, metadata_filename) 
{
  path.to.h5file <- file.path(projectpath, projectname, 'data', filename)
  path.to.metadatafile <- file.path(
    projectpath, projectname, 'data', metadata_filename)
  
  # Read out observations
  adata.obs <- read.csv(file=path.to.metadatafile)
  
  # Read out gene infos
  adata.vars <- rhdf5::h5read(path.to.h5file, "/var/")
  
  # Read out filtered, raw counts
  counts_data <- rhdf5::h5read(path.to.h5file, "/layers/counts")

  # add row and column names
  colnames(counts_data) <- adata.obs$Helmholtz_identifyer
  rownames(counts_data) <- adata.vars$hgnc_symbol
  
  rownames(adata.obs) <- adata.obs$Helmholtz_identifyer
  
  adata.vars <- convert_obs_contents_to_df(obs_list=adata.vars)
  rownames(adata.vars) <- adata.vars$hgnc_symbol
  
  
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts_data),
    colData = adata.obs, rowData = adata.vars)
  
  return(sce)
}


def_sort_proportions <- function(df)
{
  cell_order <- df %>%
    group_by(CellGroup) %>%
    summarize(avg_prop = mean(Proportion)) %>%
    arrange(avg_prop) %>%   # ascending order instead of descending
    pull(CellGroup)
  
  df_sorted <- df %>%
    mutate(CellGroup = factor(CellGroup, levels = cell_order)) %>%
    arrange(CellGroup, desc(Proportion))  # sorts within each cell type by proportion
  
  return(df_sorted)
  
}


# ============================================================================ #
#                           Set Parameters; Load data                          #
# ============================================================================ #
# Input directory
path_wd = '/Volumes/T7/Projects/Eyerich_AG_projects'
project = 'Bulk_Deconvolution__Christina_Hillig'

# Name of input files
input.h5file.sce <- 'SC_HaniffaLab_Lesion.h5' # 'SC_HaniffaLab_Lesion_with_TRM_stricter_thresholds.h5'
input.metafile.sce <- 'SC__MetaData.csv'
input.h5file.bulk <- 'Lesion_RNAseq_20210720_patient_meta_data_v04__CH__Endotypes_230620.h5'
input.metafile.bulk <- 'Bulk__MetaData.csv'

# This is a user-specific parameter which removes rare cell types
rare_ct_cut_off <- 100  

nice_colors <- c(
  "#FF0000",   # bright red
  "darkred",   # darker red
  "#D55E00",   # reddish-orange
  "darkorange",# orange
  "#E69F00",   # yellow-orange
  "#F0E442",   # bright yellow
  "yellow",    # pale yellow
  "#009E73",   # green
  "darkgreen", # darker green
  "green",     # regular green
  "#56B4E9",   # light blue
  "#0072B2",   # medium blue
  "darkblue",  # dark blue
  "purple",    # purple
  "#CC79A7",   # pink-purple
  "pink",      # bright pink
  "brown",     # brown
  "#000000",   # black
  "#999999"    # gray
)

# Output directory
save_folder = file.path(
  '/Volumes', 'CH__data', 'Projects', 'Eyerich_AG_projects',
  'BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer', 'analysis', 
  'Molecular_subtypes', 'output', 'Deconvolution_MuSiC', Sys.Date())
dir.create(save_folder, showWarnings = FALSE,  recursive = TRUE)


# Load data
adata.sce <- def_load_data(
  projectpath=path_wd, projectname=project, filename=input.h5file.sce, 
  metadata_filename=input.metafile.sce) 

adata.bulk <- def_load_bulk_data(
  projectpath=path_wd, projectname=project, filename=input.h5file.bulk, 
  metadata_filename=input.metafile.bulk) 

# # Sanity Check - subset to eczema and psoriasis
# adata.bulk <- adata.bulk[, adata.bulk$diag %in% c("psoriasis", "eczema")]

# ============================================================================ #
#                                Data preparation                              #
# ============================================================================ #
# Note Quality control by means of removing low quality cells and genes has 
# already been perform on the single-cell data. So, we skip that.

# Number of patients: 7
# Number of samples: 28 Lesion samples

# 1.1 Removing very rare cell types
# Count cells per cluster
ct_counts <- table(colData(adata.sce)$full_clustering)
# these cell types are removed: ILC2 66; Schwann_2 63; Plasma 2
ct_to_keep <- names(ct_counts[ct_counts > rare_ct_cut_off])
adata.sce <- adata.sce[, colData(adata.sce)$full_clustering %in% ct_to_keep]

# 1.2 Keep only double-positive (FOXP3 & CTLA4) Tregs 
# Identify Treg cells
treg_cells <- adata.sce$full_clustering == "Treg"
# # Extract counts for FOXP3 and CTLA4
# foxp3_counts <- assay(adata.sce, "counts")["FOXP3", ]
# ctla4_counts <- assay(adata.sce, "counts")["CTLA4", ]
# # Mask for double-positive Tregs
# doublepos_tregs <- treg_cells & (foxp3_counts > 0) & (ctla4_counts > 0)
# # Keep cells: either non-Tregs or double-positive Tregs
# keep_cells <- !treg_cells | doublepos_tregs
# # Subset the SCE object
# adata_filtered <- adata.sce[, keep_cells]

# Remove Tregs
adata_filtered <- adata.sce[, !treg_cells]
# Verify the new counts
table(adata_filtered$full_clustering)


# 2. Filter for the shared genes across bulk and single-cell data before 
# selecting highly variable genes.
shared_genes <- intersect(rownames(adata.bulk), rownames(adata_filtered))
length(shared_genes) # 16491
# Filter scRNaseq dataset to include only shared genes
adata_filtered_shared <- adata_filtered[shared_genes, ]

# 3. Normalise counts
# Estimate size factors using total counts
adata_filtered_shared <- scuttle::computeLibraryFactors(adata_filtered_shared)
# Log-normalize
adata_filtered_shared <- scuttle::logNormCounts(adata_filtered_shared)

# 4. Highly variable gene selection
hvgs <- scran::modelGeneVar(adata_filtered_shared)
# Select top 5000 highly variable genes
top_hvgs <- scran::getTopHVGs(hvgs, n = 5000)  

# Manually add genes also used in scRNAseq paper
# Create a vector of gene names
# genes <- c(
#   "ACKR1C1", "PLIN2", "BNIP3L", "FTH1", "CTSL", "FKBP1A", "EDNRB", "SQSTM1",
#   "PRKAR1A", "SOD2", "GSTO1", "CHCHD2", "ANXA5", "FTL", "PRDX1", "TUBA1C",
#   "NDUFS2", "NQO1", "CAPZB", "TPI1", "SEC61G", "ENO1", "MEDAG", "SDCBP",
#   "FAM96B", "PRMT1", "NHRNPA2B1", "PTGES", "SELENOM", "TMBIM6", "HLA-A", "MIF",
#   "PSMD8", "GYPC", "SLC43A3", "SNRPG", "NDUFC2", "HIPK2", "TUBA1B", "HINT1",
#   "PSMA7", "HMOX1", "ATP6AP2", "ATP6V0E1", "HTATIP2", "PSMB1", "TNFAIP6",
#   "WIPI1", "RGS16", "MT1A", "MYL9", "MT1M", "ZFP36", "ACTA2", "TAGLN", "S100A4",
#   "ADIRF", "SOCS3", "MYH11", "GADD45B", "HSPA1B", "HSPA1A", "HSPB1", "DUSP1",
#   "DNAJB1", "TPM1", "CRIP1", "JUNB", "ATF3", "TPM2", "FLNA", "TINAGL1", "ACTG1",
#   "BGN", "SPARC", "SH3BRGL", "RBM3", "SPARCL1", "TIMP3", "PTMS", "IGFBP7", "VIM",
#   "S100A6", "RPL37", "RPL34", "JAG1", "CHCHD10", "PPDPF", "IGFBP4", "CRTAP",
#   "AHNAK", "FRZB", "CSRP1", "FOS", "PPP1R14A", "POSTN", "RPL22", "CEBPD"
# )
genes <- c()
# # Combine HVGs + manual genes
top_hvgs_paper_genes <- union(top_hvgs, genes)
shared_hvgs <- intersect(top_hvgs_paper_genes, shared_genes)

# Subset the scRNAseq data to highly variable genes (for visualization)
adata_filtered <- adata_filtered[shared_hvgs, ]

# 5. Prepare bulk data as an ExpressionSet
# First filter bulk to also contain only shared genes and the manually added ones
adata.bulk <- adata.bulk[shared_hvgs, ]
bulk_matrix <- as.matrix(assay(adata.bulk, "counts"))
# bulk_eset <- ExpressionSet(assayData = bulk_matrix)

# ============================================================================ #
#                             Perform Deconvolution                            #
# ============================================================================ #
# Combine majory cell types
celltype_dict <- list(
  # KCs = c("Proliferating_KC", "Differentiated_KC", "Differentiated_KC*", "Undifferentiated_KC"),
  `Differentiated KCs` = c("Differentiated_KC", "Differentiated_KC*"),
  `Proliferating KCs` = c("Proliferating_KC"),
  `Undifferentiated KCs` = c("Undifferentiated_KC"),
  Fibroblasts = c("F1", "F2", "F3"),
  VEs = c("VE1", "VE2", "VE3"),
  LEs = c("LE1", "LE2"),
  Pericytes = c("Pericyte_1", "Pericyte_2"),
  Melanocytes = c("Melanocyte"),
  `Mast cells` = c("Mast_cell"),
  `Mø` = c("Macro_1", "Macro_2", "Inf_mac"),
  `APC` = c("DC1", "DC2", "LC_1", "LC_2", "LC_3", "LC_4",
            "MigDC", "Mono_mac", "moDC_1", "moDC_2", "moDC_3"),
  `NK/ILC` = c("NK", "ILC1_NK", "ILC1_3", "ILC2"),
  `Schwann cells` = c("Schwann_1", "Schwann_2"),
  `T cells` = c("Tc", "Tc_IL13_IL22", "Tc17_Th17", "Th", "Treg", "T_RM"),
  # `Tc` = c("Tc"),
  # `Tc_IL13_IL22` =c("Tc_IL13_IL22"),
  # `Tc17_Th17` = c("Tc17_Th17"),
  # `Treg` = c("Treg"),
  # `Th` = c("Th"),
  `Plasma cells` = c('Plasma')
)

# ---> Aggregate to Major Cell types
# # Flatten mapping
# cluster_to_major <- unlist(lapply(names(celltype_dict), function(x) {
#   setNames(rep(x, length(celltype_dict[[x]])), celltype_dict[[x]])
# }))
# 
# # Map each cell to major type
# adata.sce$major_type <- cluster_to_major[adata.sce$full_clustering]
# table(adata.sce$major_type)

# # ----> Merge only T-cells
# # Create a new cluster column for deconvolution
# colData(adata.sce)$cluster_for_deconv <- colData(adata.sce)$full_clustering
# # Convert factor to character
# colData(adata.sce)$cluster_for_deconv <- as.character(colData(adata.sce)$cluster_for_deconv)
# # Combine T_RMs with the main T-cell cluster
# tcell_clusters <- c("Tc", "Tc_IL13_IL22", "Tc17_Th17", "Th", "Treg", "T_RM")
# colData(adata.sce)$cluster_for_deconv[
#   colData(adata.sce)$full_clustering %in% tcell_clusters
# ] <- "T cells"
# colData(adata.sce)$cluster_for_deconv <- factor(colData(adata.sce)$cluster_for_deconv)

# Use only raw counts -> remove the logcounts assay
assay(adata_filtered, "logcounts") <- NULL

# add patient ID
colData(adata_filtered)$patient_id <- paste0(
  adata_filtered$Status, "_", adata_filtered$Age
)

# Estimate cell type proportions
# This function is to calculate the MuSiC deconvolution proportions
est_prop.bulk <- MuSiC::music_prop(
  bulk.mtx = bulk_matrix,
  sc.sce = adata_filtered,
  clusters = 'full_clustering', 
  samples = 'patient_id', 
  # select.ct = c("Tc", "Tc_IL13_IL22", "Tc17_Th17", "Th", "Treg", 
  #               "NK", "ILC1_NK", "ILC1_3", "ILC2"),
  verbose = F)
# Save result
saveRDS(est_prop.bulk, file.path(
  save_folder, "Estimated_Celltype_probabilities.rds"))

# View the estimated proportions
head(est_prop.bulk$Est.prop.weighted)

# Step 1: Extract proportionsm and convert cell type proportions to a data frame with sample names
props <- est_prop.bulk$Est.prop.weighted 
# # Fix the NA column name
# colnames(props)[is.na(colnames(props))] <- "Unknown"

props <- props %>%
  as.data.frame() %>%
  rownames_to_column(var = "Sample")

# Select columns of interest in bulk metadata
bulk_meta <- as.data.frame(colData(adata.bulk))
bulk_meta <- bulk_meta[c('Endotypes', 'Helmholtz_identifyer', "diag")]

# Step 2: Merge with your sample metadata (contains Endotype info)
all(props$Sample == bulk_meta$Helmholtz_identifyer)
props <- left_join(props, bulk_meta, by = c("Sample" = "Helmholtz_identifyer"))


# Create a collapsed data.frame, keeping only the sample ID and Endotype columns
props_collapsed <- props[, c("Sample", "Endotypes", "diag")]


# Add one column per grouped cell type
for (group in names(celltype_dict)) {
  celltypes <- intersect(celltype_dict[[group]], colnames(props))
  props_collapsed[[group]] <- rowSums(props[, celltypes, drop = FALSE])
}
# =====================

# save to xlsx or csv
openxlsx::write.xlsx(x=props_collapsed, 
                     file=file.path(save_folder, 'Celltype_deconvolution.xlsx'))


# ============================================================================ #
#                    Plot cell type composition in endotypes                   #
# ============================================================================ #
# Average per cell group per endotype
props_endotype <- props_collapsed  %>% # props props_collapsed
  select(-Sample, -diag) %>%
  group_by(Endotypes) %>%
  summarise(across(everything(), mean))

# Convert proportion matrix to long format
props_long <- props_endotype %>%
  pivot_longer(-Endotypes, names_to = "CellGroup", values_to = "Proportion")

props_long$Endotypes <- factor(
  props_long$Endotypes, levels = c(
    "E5", "E6", "E7", "E8", "E9", "E10",
    "E11", "E12", "E13",
    "E1", "E2", "E3", "E4"))

# Remove those cell types with 0 probabilities
props_long$CellGroup <- as.factor(props_long$CellGroup)
props_long <- props_long %>%
  mutate(CellGroup=forcats::fct_reorder(CellGroup, Proportion),.desc=T)

# Remove cell types with 0 Proportion
props_long <- props_long %>% filter(Proportion > 0)

props_long_sorted <- def_sort_proportions(props_long)

# Step 5: Plot stacked barplot of proportions grouped by Endotype
p.barplot.celltypes <- ggplot(
  props_long_sorted, aes(x = Endotypes, y = Proportion, fill = CellGroup)) +
  geom_bar(stat = "identity", position = "fill") +
  # scale_fill_brewer(palette = "Set3")
  scale_fill_manual(values = nice_colors, 
    name='Cell types') + 
  theme_minimal() +
  ylab("Proportion") +
  labs(title = "Cell type composition per Endotype")

pdf(file = file.path(save_folder, paste0("Barplot_celltypes_endotypes.pdf")), 
    width = 6, height=4)
print(p.barplot.celltypes)
dev.off()


# ============================================================================ #
#                   Plot cell type composition in diagnosis                    #
# ============================================================================ #


props_diag <- props_collapsed  %>% # props props_collapsed
  select(-Sample, -Endotypes) %>%
  group_by(diag) %>%
  summarise(across(everything(), mean))

props_long_diag <- props_diag %>%
  pivot_longer(-diag, names_to = "CellGroup", values_to = "Proportion")

# Remove cell types with 0 Proportion
props_long_diag <- props_long_diag %>% filter(Proportion > 0)

props_long_diag_sorted <- def_sort_proportions(props_long_diag)


p.celltypes_diag <- ggplot(
  props_long_diag_sorted, aes(x = diag, y = Proportion, fill = CellGroup)) +
  geom_bar(stat = "identity", position = "fill") +
  # scale_fill_brewer(palette = "Set3")
  scale_fill_manual(values = nice_colors, 
    name='Cell types') + 
  theme_minimal() +
  ylab("Proportion") +   xlab("") +
  labs(title = "Cell type composition per disease") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  guides(fill = guide_legend(ncol = 2, title.position = "top"))  # Split legend into 2 columns

pdf(file = file.path(save_folder, paste0("Barplot_celltypes_diagnosis.pdf")), 
    width = 10, height=6)
print(p.celltypes_diag)
dev.off()


# ============================================================================ #
#                       Plot immune cell composition                           #
# ============================================================================ #
pastel_colors <- c(
  "#FFB3BA",  # soft red/pink
  "#FFDFBA",  # soft orange
  "#FFFFBA",  # soft yellow
  "#BAFFC9",  # soft green
  "#BAE1FF",  # soft blue
  "#E0BAFF",  # soft purple/lavender
  "#FED9A6", # pastel orange
  "#CFCFCF",  # pastel gray
  "#FFC0CB",  # pink
  "#FFDAB9",  # peach
  "#FFFACD",  # lemon
  "#C1E1C1",  # mint
  "#ADD8E6",  # light blue
  "#D8BFD8"   # thistle
)


# show T-cell composition: No 'ILC2'
required_cols <- c(
  "Sample", "Endotypes", "diag", "Tc", "Tc_IL13_IL22", "Tc17_Th17", "Th", "Treg", 
  'ILC2', "NK", "ILC1_NK", "ILC1_3", "T_RM")
# Keep only those columns that exist in props
props_tcells <- props[, intersect(required_cols, colnames(props))]


# Average per cell group per endotype
tcells_props_endotype <- props_tcells %>%
  select(-Sample, -diag) %>%
  group_by(Endotypes) %>%
  summarise(across(everything(), mean))

tcells_props_long <- tcells_props_endotype %>%
  pivot_longer(-Endotypes, names_to = "CellGroup", values_to = "Proportion")

tcells_props_long$Endotypes <- factor(
  tcells_props_long$Endotypes, levels = c(
    "E5", "E6", "E7", "E8", "E9", "E10",
    "E11", "E12", "E13",
    "E1", "E2", "E3", "E4"))

tcells_props_long$CellGroup <- as.factor(tcells_props_long$CellGroup)
tcells_props_long <- tcells_props_long %>%
  mutate(CellGroup=forcats::fct_reorder(CellGroup, Proportion),.desc=T) 

# Remove cell types with 0 Proportion
tcells_props_long <- tcells_props_long %>% filter(Proportion > 0)

tcells_props_long_sorted <- def_sort_proportions(tcells_props_long)

p.barplot.tcells <- ggplot(
  tcells_props_long_sorted, aes(x = Endotypes, y = Proportion, fill = CellGroup)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_brewer(palette = "Set3", name='Cell types') +
  scale_fill_manual(values = pastel_colors)+
  theme_minimal() + 
  ylab("Proportion") +
  labs(title = "Immune cell type composition per Endotype", 
       fill = "Immune Cell Type")

pdf(file = file.path(save_folder, paste0("Barplot_Immune_cells_endotypes.pdf")), 
    width = 6, height=4)
print(p.barplot.tcells)
dev.off()


# Average per cell group per diagnosis
tcells_props_diag <- props_tcells %>%
  select(-Sample, -Endotypes) %>%
  group_by(diag) %>%
  summarise(across(everything(), mean))

tcells_props_diag_long <- tcells_props_diag %>%
  pivot_longer(-diag, names_to = "CellGroup", values_to = "Proportion")

tcells_props_diag_long$CellGroup <- as.factor(tcells_props_diag_long$CellGroup)
tcells_props_diag_long <- tcells_props_diag_long %>%
  mutate(CellGroup=forcats::fct_reorder(CellGroup, Proportion),.desc=T) 

# Remove cell types with 0 Proportion
tcells_props_diag_long <- tcells_props_diag_long %>% filter(Proportion > 0)

tcells_props_diag_long_sorted <- def_sort_proportions(tcells_props_diag_long)

p.barplot.tcells_diag <- ggplot(
  tcells_props_diag_long_sorted, aes(x = diag, y = Proportion, fill = CellGroup)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_brewer(palette = "Set3", name='Cell types') +
  scale_fill_manual(values = pastel_colors)+
  theme_minimal() + 
  ylab("Proportion") + xlab("") +
  labs(title = "Immune cell type composition per diagnosis", 
       fill = "Immune Cell Type") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf(file = file.path(save_folder, paste0("Barplot_Immune_cells_diagnosis.pdf")), 
    width = 6, height=6)
print(p.barplot.tcells_diag)
dev.off()


# # ============================================================================ #
# #            Estimate cell-type-specific expression per sample                 #
# # ============================================================================ #
# # Need tools like CIBERSORTx HiRes
# 
# bulk_matrix <- est_prop.bulk$Est.prop.weighted
# tcell_subtypes <- c("Tc", "Tc17_Th17", "Tc_IL13_IL22", "Th", "Treg", "T_RM")
# 
# # Initialize output with same rownames as bulk_matrix
# bulk_tcell_subtypes <- data.frame(matrix(nrow = nrow(bulk_matrix), ncol = 0))
# rownames(bulk_tcell_subtypes) <- rownames(bulk_matrix)
# 
# for (subtype in tcell_subtypes) {
#   fraction_within_T <- sum(colData(adata.sce)$full_clustering == subtype) /
#     sum(colData(adata.sce)$cluster_for_deconv == "T cells")
#   
#   # Ensure we are multiplying the correct column
#   bulk_tcell_subtypes[, subtype] <- bulk_matrix[, "T cells", drop = TRUE] * fraction_within_T
# }
# 
# # Make bulk_tcell_subtypes a data.frame with Sample column
# bulk_tcell_subtypes$Sample <- rownames(bulk_tcell_subtypes)
# 
# # Merge with bulk_meta to get Endotypes
# bulk_tcell_subtypes <- left_join(
#   bulk_tcell_subtypes,
#   bulk_meta,
#   by = c("Sample" = "Helmholtz_identifyer")
# )
# 
# 
# 
# bulk_tcell_subtypes <- est_prop.bulk$Est.prop.weighted[, "T cells" , drop = FALSE] %>% 
#   as.data.frame() %>% 
#   rownames_to_column("Sample")
# 
# for(subtype in tcell_subtypes){
#   # Fraction within T cells in the scRNA-seq reference
#   fraction_within_T <- sum(colData(adata.sce)$full_clustering == subtype) / 
#     sum(colData(adata.sce)$cluster_for_deconv == "T cells")
#   
#   # Multiply per sample
#   bulk_tcell_subtypes[[subtype]] <- bulk_tcell_subtypes$`T cells` * fraction_within_T
# }
# 
# # Add Endotype info
# bulk_meta <- as.data.frame(colData(adata.bulk))[ , c("Endotypes", "Helmholtz_identifyer")]
# bulk_tcell_subtypes <- merge(bulk_tcell_subtypes, bulk_meta,
#                              by.x = "Sample", by.y = "Helmholtz_identifyer",
#                              all.x = TRUE)
# 
# # Convert to long format for plotting
# bulk_tcell_long <- bulk_tcell_subtypes %>%
#   pivot_longer(
#     cols = all_of(tcell_subtypes), 
#     names_to = "TcellSubtype", 
#     values_to = "Proportion"
#   )
# 
# p.barplot.tcells <- ggplot(
#   bulk_tcell_long, aes(x = Endotypes, y = Proportion, fill = TcellSubtype)) +
#   geom_bar(stat = "identity", position = "fill") +
#   scale_fill_brewer(palette = "Set3", name='Cell types')
# 
# 
# tcell_fractions_per_sample <- colData(adata.sce) %>%
#   as.data.frame() %>%
#   filter(cluster_for_deconv == "T cells") %>%       # keep only T cells
#   group_by(sample_id, full_clustering) %>%
#   summarise(n = n(), .groups = 'drop') %>%
#   group_by(sample_id) %>%
#   mutate(Fraction = n / sum(n)) %>%                # fraction within T cells per sample
#   filter(full_clustering %in% tcell_subtypes)     # keep only the T-cell subtypes you care about
# 
# 
# # bulk_tcells: dataframe with columns Sample, T_cells (total T-cell proportion)
# # Extract bulk T-cell proportions per sample
# bulk_tcells <- est_prop.bulk$Est.prop.weighted[ , "T cells", drop = FALSE] %>%
#   as.data.frame() %>%
#   rownames_to_column("Sample") %>%
#   rename(T_cells = "T cells")
# 
# bulk_tcell_subtypes <- tcell_fractions_per_sample %>%
#   left_join(bulk_tcells, by = c("sample_id" = "Sample")) %>%
#   mutate(AbsoluteProportion = T_cells * Fraction) %>%
#   rename(TcellSubtype = full_clustering, Sample = sample_id)
# 
# bulk_tcell_subtypes <- bulk_tcell_subtypes %>%
#   left_join(bulk_meta, by = c("Sample" = "Helmholtz_identifyer"))
# 
# # Already long if using above, otherwise:
# bulk_tcell_long <- bulk_tcell_subtypes %>%
#   select(Sample, Endotypes, TcellSubtype, AbsoluteProportion)
# 
# ggplot(
#   bulk_tcell_long, aes(x = Endotypes, y = AbsoluteProportion, fill = TcellSubtype)) +
#   geom_bar(stat = "identity", position = "fill") +
#   scale_fill_brewer(palette = "Set3", name='Cell types')





