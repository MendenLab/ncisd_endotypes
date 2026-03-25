# Install packages

renv::init()

proj_lib <- file.path(getwd(), "r_libs")
dir.create(proj_lib, showWarnings = FALSE)
# Add to library path
.libPaths("r_libs")

# check permissions
file.access(proj_lib, 2) == 0  # TRUE means writable


# CRAN packages
pkgs <- c(
  "BiocManager", "openxlsx", "forcats", "tibble", "ggplot2",
  "data.table", "stringr", "dplyr", "plyr", "cowplot", 'tidyverse'
)
install.packages(pkgs, lib = proj_lib)


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager", lib = proj_lib)

# Order is important
pkgs <- c("rhdf5", "scuttle", "scater", "scran", "SingleCellExperiment")

for (p in pkgs) {
  if (p %in% rownames(installed.packages(lib.loc = proj_lib))) {
    remove.packages(p, lib = proj_lib)
  }
}

BiocManager::install(
  c("rhdf5", "scuttle", "scater", "scran", "SingleCellExperiment"),
  lib = proj_lib,
  force = TRUE
)
# Verify installation
sapply(c("rhdf5", "scuttle", "scater", "scran", "SingleCellExperiment"),
       function(p) system.file(package = p, lib.loc = proj_lib))


# install.packages('devtools', lib = proj_lib)
# # Had to reinstall rlang, lifecycle, glue, cli -> better skip updates
# library(devtools)
# devtools::install_github('xuranw/MuSiC', lib = proj_lib)
# devtools::install_github("meichendong/SCDC", lib = proj_lib)
renv::install('xuranw/MuSiC', lib = proj_lib)


renv::snapshot()

