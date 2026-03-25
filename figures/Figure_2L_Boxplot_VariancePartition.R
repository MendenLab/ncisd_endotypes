rm(list = ls(all=TRUE))
library(variancePartition)
library(lsr)

######################################################################
#                                Start                               #
######################################################################
data.root <- file.path('/Volumes/CH__data/Projects/Eyerich_AG_projects', 
                       'BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer')
input.dir <- file.path(data.root, 'analysis/Molecular_subtypes', 
                        'output/Variance_explained', 
                        'variationPartition/variancePartition_input', 
                       '2023-11-23')
output.dir <- file.path(data.root, 'analysis/Molecular_subtypes', 
                       'output/Figure_3K_Varianzpartition', Sys.Date())
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

vp_batch_corrected <- readRDS(file = file.path(input.dir, 'varPart_Endotypes_batchcorrected.rds'))

# Plot preparation
df.reduced <- vp_batch_corrected
df.reduced$Residuals <- NULL
data.table::setnames(df.reduced, "diag", "Diagnosis")
data.table::setnames(df.reduced, "Molecular.subtype.res0.9", "Endotypes")

df.reduced.melted <- reshape2::melt(df.reduced)

plotPercentBars( df.reduced[1:10,] )
plotVarPart( df.reduced , col=c('aliceblue', 'aliceblue', 'aliceblue'), 
             label.angle = 90)

df.reduced.melted$variable <- as.factor(df.reduced.melted$variable)
df.reduced.melted$variable <- factor(
  df.reduced.melted$variable, levels = c('Diagnosis', 'Pattern', 'Endotypes'))
df.reduced.melted$value <- df.reduced.melted$value * 100

# Plot
pdf(file = file.path(output.dir,  
                     paste0('Variance_explained_batch_corrected', '.pdf')), 
    width = 6, height = 4)
p <- ggplot(df.reduced.melted, aes(x=variable, y=value, fill=variable)) + 
  geom_boxplot() + labs(x="", y = "Variance explained (%)")
p <- p + scale_fill_grey() + theme_classic()  + 
  theme(legend.position = "none", text=element_text(size = 20))
print(p)
dev.off()


# Calculate mean and median 
mean(as.numeric(df.reduced$Endotypes) * 100) # 11.54916
mean(as.numeric(df.reduced$Diagnosis) * 100) # 4.953487
mean(as.numeric(df.reduced$Pattern) * 100) # 1.851041

median(as.numeric(df.reduced$Endotypes) * 100) # 8.699944
median(as.numeric(df.reduced$Diagnosis) * 100) # 1.386536
median(as.numeric(df.reduced$Pattern) * 100) # 0

sd(as.numeric(df.reduced$Endotypes) * 100) # 10.57167
sd(as.numeric(df.reduced$Diagnosis) * 100) # 7.676689
sd(as.numeric(df.reduced$Pattern) * 100) # 4.426373

# Calculate effect size - should not be used when sample are identical in each group
cohen_d_ED <- lsr::cohensD(df.reduced$Endotypes * 100 , df.reduced$Diagnosis * 100) #  0.7139507 (medium effect)
cohen_d_EP <- lsr::cohensD(df.reduced$Endotypes * 100, df.reduced$Pattern * 100) # 1.196693 (strong effect)
cohen_d_PD <- lsr::cohensD(df.reduced$Pattern * 100, df.reduced$Diagnosis * 100) # 0.4951275 (low effect)




