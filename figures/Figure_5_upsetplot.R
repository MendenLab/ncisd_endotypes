set.seed(1)

library(stringr)
library(UpSetR)
library(ComplexHeatmap)
library(xlsx)

classifier <- 'IL23_IL17' # 'IL23_IL17'  # 'TNF'

if (classifier == 'TNF') 
{
  # df <- read.csv('/Volumes/CH__data/Projects/Eyerich_AG_projects/BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer/analysis/Molecular_subtypes/output/Compare_GEx_TNF_models/2024-01-09/UpSetPlot_features_modelnames.csv')
  df <- read.csv('/Volumes/CH__data/Projects/Eyerich_AG_projects/BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer/analysis/Molecular_subtypes/output/Compare_GEx_TNF_models/2024-05-06/UpSetPlot_features_modelnames.csv')
  
} else 
{
  # df <- read.csv('/Volumes/CH__data/Projects/Eyerich_AG_projects/BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer/analysis/Molecular_subtypes/output/Compare_GEx_IL23_IL17_models/2023-12-28/UpSetPlot_features_modelnames.csv')
  df <- read.csv('/Volumes/CH__data/Projects/Eyerich_AG_projects/BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer/analysis/Molecular_subtypes/output/Compare_GEx_IL23_IL17_models/2024-05-06/UpSetPlot_features_modelnames.csv')
}

features <- gsub("\\[|\\]", "", df$features)
features <- gsub("'", "", features)
named_list.features <- str_split(features, ', ')

names(named_list.features) <- df$model


# intersection.features <- UpSetR::fromList(named_list.features)
# 
# UpSetR::upset(data=intersection.features, order.by = "freq", 
#               nsets = dim(df)[1], keep.order=T, sets = rev(as.vector(df$model)),
#               nintersects = NA)



savefolder <- file.path('/Volumes/CH__data/Projects/Eyerich_AG_projects/BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer/analysis/Molecular_subtypes/output/Figure_5', 
                        Sys.Date())
dir.create(savefolder, recursive = TRUE, showWarnings = FALSE)

m <- ComplexHeatmap::list_to_matrix(named_list.features)
xlsx::write.xlsx(m, file.path(savefolder, paste0('Figure_5_', classifier, '_sets.xlsx')))
m = ComplexHeatmap::make_comb_mat(m, mode = 'distinct')


pdf(file = file.path(savefolder, 
                     paste0('Figure_5_', classifier, '_Upsetplot.pdf')), 
    width = 6, height = 6)
p <- ComplexHeatmap::UpSet(
  t(m), set_order =  as.vector(df$model), pt_size = unit(2, "mm"), lwd = 3)
  # comb_col = c("red", "grey"))
print(p)
dev.off()


# For figure 
combination_names <- comb_name(m)
set_size(m)
comb_size(m)

# Initialize an empty list to hold columns of different lengths
columns_list <- list()
for (sets in combination_names) 
{
  # Add the column data to the list
  # returns features shared across combination set sets
  columns_list[[sets]] <- extract_comb(m, sets)
}

# Merge columns of different lengths into a data frame
df_sets_fs <- as.data.frame(lapply(columns_list, `length<-`, max(lengths(columns_list))))
# save feature combination sets
xlsx::write.xlsx(
  t(df_sets_fs),  file.path(savefolder, 
                         paste0('Figure_5_', classifier, '_feature_sets.xlsx')))

