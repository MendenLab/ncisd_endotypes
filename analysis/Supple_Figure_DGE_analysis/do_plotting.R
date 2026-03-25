rm(list = ls(all=TRUE))

library('xlsx')
library('magrittr')
library(data.table)

general_path = '/Volumes/CH__data/Projects/Eyerich_AG_projects/BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer/analysis/Molecular_subtypes'

source(file.path(general_path, 'r_scripts', 'DGE_analysis', 'Volcanoplot.R'))

# 1. Create output directory
save_dir <- file.path(general_path, 'output', 'Figure_S3_Volcanoplots', Sys.Date())
dir.create(save_dir, showWarnings = FALSE,  recursive = TRUE)

# 2. Set .xlsx files
input.dir.degs <- file.path(general_path, 'output','DGE_analysis', 
                            'Endotypes_Dendrogram', '2023-08-29')

file = 'DEGs_E1_E2__ vs __E3_E4__Age_Sex_batchID_subtypes.xlsx'

group_one <- strsplit(strsplit(file, '__')[[1]][1], "s_")[[1]][2]
group_two <- strsplit(file, '__')[[1]][3]

comparison.name <- paste0(group_one, '_vs_', group_two)

# Load DGE list
df_degs <- xlsx::read.xlsx(file.path(input.dir.degs, file), 
  sheetIndex = 1, header = TRUE, 
  Comparision=comparison.name, )

# Add labels to E1/E2 vs E3/E4
labels = c('IL6', 'IL1B')

# 5. Plot Volcano with marker genes
pdf(file = file.path(save_dir, paste0(comparison.name, '_Volcanoplot', '.pdf')), 
    width = 6, height = 6)
print(
  def.get.volcanoplot(df.data=df_degs, log2fc.cut=1, add.labels=labels) 
)
dev.off()