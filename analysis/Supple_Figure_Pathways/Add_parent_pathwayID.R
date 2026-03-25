library('xlsx')

df.test <- read.csv('/Users/christina.hillig/Downloads/Ensembl2Reactome_All_Levels.txt', sep='\t')
df.relation <- read.csv('/Users/christina.hillig/Downloads/ReactomePathwaysRelation-2.txt', sep='\t')
# df.human_hierarchy <- fread('/Users/christina.hillig/Downloads/Complex_2_Pathway_human.txt')

colnames(df.test) <- c("ENSEMBLID", "ID", "Weblink", "Description", "Level", "Species")
colnames(df.relation) <- c("ParentPathway", "ChildPathway")

df.test <- df.test[df.test$Species == 'Homo sapiens', ]


general_path = '/Volumes/CH__data/Projects/Eyerich_AG_projects/BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer/analysis/Molecular_subtypes'
general.dir <- file.path(general_path, 'output/Pathways/Endotypes_Dendrogram/2024-12-10')

# Use list.files to find all files in subdirectories
all_files <- list.files(general.dir, recursive = TRUE, full.names = TRUE, pattern = "\\.xlsx")


# Identify parent pathways in Dendrogram splits
all.parent.descriptions <- c()
for (file_name in all_files) {
  df.pas <- read.xlsx(file_name, sheetIndex = 1)
  
  if (nrow(df.pas) > 0) 
  {
    pa.ids <- df.pas$ID
    parent.descriptions <- c()
    for (id in pa.ids) {
      # Get parent pathway name
      parent.id <- df.relation[which(df.relation$ChildPathway %in% id), 'ParentPathway']
      parent.description <- df.test[which(df.test$ID %in% parent.id), 'Description'][1]
      parent.descriptions <- c(parent.descriptions, parent.description)
    }
    df.pas$parentID <- parent.descriptions
    write.xlsx(x=df.pas, file=file_name)
    
    all.parent.descriptions <- c(all.parent.descriptions, parent.descriptions) 
  }
  
}



# ============ Identify parent pathways in 1 vs Rest ============ #
general_path = '/Volumes/CH__data/Projects/Eyerich_AG_projects/BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer/analysis/Molecular_subtypes'
general.dir.1vsRest <- file.path(general_path, 'output/Pathways/Endotypes_1_vs_R_ORA_GSEA/2024-12-10')
# general.dir.1vsRest.ora <- file.path(general_path, 'output/Pathways/Molecular_subtypes_res0.9_1vsRest/unique_Pas/2024-07-16')

# Use list.files to find all files in subdirectories
all_files.gsea <- list.files(general.dir.1vsRest, recursive = TRUE, full.names = TRUE, pattern = "\\.xlsx")
# all_files.ora <- list.files(general.dir.1vsRest.ora, recursive = TRUE, full.names = TRUE)


all.parent.descriptions.1vsRest.gsea <- c()
for (file_name in all_files.gsea) {
  df.pas <- read.xlsx(file_name, sheetIndex = 1)
  
  if (nrow(df.pas) > 0) 
  {
    pa.ids <- df.pas$ID
    parent.descriptions <- c()
    for (id in pa.ids) {
      # Get parent pathway name
      parent.id <- df.relation[which(df.relation$ChildPathway %in% id), 'ParentPathway']
      parent.description <- df.test[which(df.test$ID %in% parent.id), 'Description'][1]
      parent.descriptions <- c(parent.descriptions, parent.description)
    }
    df.pas$parentID <- parent.descriptions
    write.xlsx(x=df.pas, file=file_name)
    
    all.parent.descriptions.1vsRest.gsea <- c(
      all.parent.descriptions.1vsRest.gsea, parent.descriptions) 
  }
  
}

# 
# all.parent.descriptions.1vsRest.ora <- c()
# for (file_name in all_files.ora) {
#   df.pas <- read.xlsx(file_name)
#   
#   if (nrow(df.pas) > 0) 
#   {
#     pa.ids <- df.pas$ID
#     parent.descriptions <- c()
#     for (id in pa.ids) {
#       # Get parent pathway name
#       parent.id <- df.relation[which(df.relation$ChildPathway %in% id), 'ParentPathway']
#       parent.description <- df.test[which(df.test$ID %in% parent.id), 'Description'][1]
#       parent.descriptions <- c(parent.descriptions, parent.description)
#     }
#     df.pas$parentID <- parent.descriptions
#     write.xlsx(x=df.pas, file=file_name)
#     
#     all.parent.descriptions.1vsRest.ora <- c(
#       all.parent.descriptions.1vsRest.ora, parent.descriptions) 
#   }
#   
# }


# Count how often a pathway pops up on the left and right side of the Dendrogram

