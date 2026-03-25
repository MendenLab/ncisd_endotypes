library(ggplot2)
library(tidyverse)
library(tidyr)
library(reshape2)


######################################################################
#                                Start                               #
######################################################################
data.root <- file.path('/Volumes/CH__data/Projects/Eyerich_AG_projects', 
                       'BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer')
# df <- readRDS(file='/Volumes/CH__data/Projects/Eyerich_AG_projects/BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer/results/variancePartition_04.04.23/varPartResid.res0.9.rds')
df <- readRDS(file.path(data.root, 'coworkers/Gülce/variancePartition_04.04.23/varPartResid.res0.9.rds'))

output.dir <- file.path(data.root, 'analysis', 'Molecular_subtypes', 'output')
output.dir <- file.path(output.dir, 'Variance_explained', 'Figure_2D', Sys.Date())
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)


df.reduced <- df
df.reduced$Residuals <- NULL
data.table::setnames(df.reduced, "diag", "Diagnosis")
data.table::setnames(df.reduced, "Molecular.subtype.res0.9", "Endotypes")

df.reduced.melted <- reshape2::melt(df.reduced)

plotPercentBars( df[1:10,] )
plotVarPart( df.reduced , col=c('aliceblue', 'aliceblue', 'aliceblue'), 
             label.angle = 90)

# outlier.colour="black", outlier.shape=16, outlier.size=1, notch=FALSE

pdf(file = file.path(output.dir,  paste0('Variance_explained', '.pdf')), 
    width = 6, height = 4)
p <- ggplot(df.reduced.melted, aes(x=variable, y=value, fill=variable)) + 
  geom_boxplot() + labs(x="", y = "Variance explained (%)")
p <- p + scale_fill_grey() + theme_classic()  + 
  theme(legend.position = "none", text=element_text(size = 20))
print(p)
dev.off()
