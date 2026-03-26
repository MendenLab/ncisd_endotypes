# Non-communicable inflammatory skin diseases comprise clinically meaningful distinct endotypes as identified by non-hypothesized integration of phenotypic and transcriptome data

Natalie Garzorz-Stark†, Christina Hillig†, Heydar Maboudi Afkham†, Martin Meinel, 
Manja Jargosch, Peter Seiringer, Felix Lauffer, Anna C Pilz, Xixi Li, Jigyansa Mishra, 
Matthias Hübenthal, Stephan Weidinger, Benjamin Klein, Lam C Tsoi, Johann Gudjonsson, 
Curdin Conrad12, Michel Gilliet, Fabian Theis, Carsten B. Schmidt-Weber, Tilo Biedermann, 
Michael P. Menden*, †, Stefanie Eyerich†, Kilian Eyerich*


*Key words*: Inflammatory skin diseases (ncISD), endotypes, precision medicine, 
disease heterogeneity, IL-23 signaling, targeted therapy, psoriasis 

## Environment
In order to run the analysis you will have to create a conda 
environment from the py38_molecularsubtypes_20230411.yml file. 
Before doing so, you will have to manually set the prefix variable in the *.yml file to your directory. 
Now, you can run the following commands in the terminal:
```bash
# The following command will create an env with the name py38_molecularsubtypes
conda env create -f py38_molecularsubtypes_20230411.yml

# Activate the conda env with
conda activate py38_molecularsubtypes_20230411
```

For the deconvolution using MuSiC the R version 4.5.1 (2025-06-13) -- "Great Square Root" is required. 
Details on which packages are needed to run the deconvolution can be found in analysis/MuSiC_deconvolution/install_pacakges.R

For the Variance Partition Analysis, DGE analysis, and pathway enrichment analysis the R version R 4.2.1 was used.
```R
# cran packages
rlibs <- c("dplyr", "gtools", "hash", "kableExtra", "knitr", "stringr", "tibble", "xlsx", "dendextend")
invisible(lapply(rlibs, require, character.only = TRUE))

# Bioconductor packages (Bioconductor version 3.12 (BiocManager 1.30.12))
bioclibs <- c("variancePartition", "DESeq2", "limma", "edgeR", "pathview",  "org.Hs.eg.db", "ReactomePA",  "enrichplot", "clusterProfiler", "DOSE", "UpSetR", "fgsea")
invisible(lapply(bioclibs, require, character.only = TRUE))
```


## Analysis

### Cluster-Informative Feature Genes (CIFG)

### Discriminative Feature Genes (DFG)

### Endotype classifier
