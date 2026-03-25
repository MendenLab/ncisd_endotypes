rm(list = ls(all=TRUE))
# Set seed before loading packages 
set.seed(1)

# Load packages
library('DESeq2')
# library(fastICA)
library(dplyr)
library(ggplot2)
library("gridExtra")
# library(factoextra)
# library(glmpca)
library(xlsx)

# library(ggsignif)
library(ggrepel)

general.dir <- '/Volumes/CH__data/Projects/Eyerich_AG_projects/BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer'


df_michigan <- read.xlsx(file.path(
  general.dir, 'analysis/Molecular_subtypes/output/DGE_analysis/Michigan_L_vs_NL/2024-01-14/DEGs_psoriasis__Age_Sex_sampleType.xlsx'), 
  sheetIndex = 1)
  
df_eyerich <- read.xlsx(file.path(
  general.dir, 'analysis/Molecular_subtypes/output/DGE_analysis/Psoriasis_L_vs_NL/2024-01-13/DEGs_psoriasis__Age_Sex_batchID_sampleType.xlsx'),
  sheetIndex = 1)
  
  
output.dir <- file.path(general.dir, 'analysis', 'Molecular_subtypes', 
                        'output', 'DGE_analysis', 'Michigan_vs_Eyerich', Sys.Date())
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
  
  
  
# --------------- signed padj plot
# create signed p-value plots
axis.title <- 16
axis.text <- 14
legend.title <- 14
legend.text <- 12

df_michigan$sign.padj.pso <- -sign(df_michigan$log2FoldChange) * log10(df_michigan$padj)
df_eyerich$sign.padj.ae <- -sign(df_eyerich$log2FoldChange) * log10(df_eyerich$padj)


df.Michigan_vs_Eyerich <- merge(
  data.frame(df_michigan[ , c('hgnc_symbol', 'sign.padj.pso')], row.names=NULL),
  data.frame(df_eyerich[ , c('hgnc_symbol', 'sign.padj.ae')], row.names=NULL),
  by = 'hgnc_symbol', all = TRUE)
df.Michigan_vs_Eyerich <- df.Michigan_vs_Eyerich[df.Michigan_vs_Eyerich$hgnc_symbol != "", ]


mask <- df.Michigan_vs_Eyerich$sign.padj.ae > -log10(0.05) & df.Michigan_vs_Eyerich$sign.padj.pso > -log10(0.05)
df.Michigan_vs_Eyerich$diffexpr <- 'Not sig.'
df.Michigan_vs_Eyerich$diffexpr[mask] <- 'Up'
mask <- df.Michigan_vs_Eyerich$sign.padj.ae < -log10(0.05) & df.Michigan_vs_Eyerich$sign.padj.pso < -log10(0.05)
df.Michigan_vs_Eyerich$diffexpr[mask] <- 'Down'
mask <- df.Michigan_vs_Eyerich$sign.padj.ae < -log10(0.05) & df.Michigan_vs_Eyerich$sign.padj.pso > -log10(0.05)
df.Michigan_vs_Eyerich$diffexpr[mask] <- 'Up Michigan - Down Eyerich'
mask <- df.Michigan_vs_Eyerich$sign.padj.ae > -log10(0.05) & df.Michigan_vs_Eyerich$sign.padj.pso < -log10(0.05)
df.Michigan_vs_Eyerich$diffexpr[mask] <- 'Up Eyerich - Down Michigan'

df.Michigan_vs_Eyerich$diffexpr <- as.factor(df.Michigan_vs_Eyerich$diffexpr)

res.Michigan_vs_Eyerich <- cor.test(
  df.Michigan_vs_Eyerich$sign.padj.pso, df.Michigan_vs_Eyerich$sign.padj.ae,
  method = "spearman", use = "complete.obs")



pdf(file = file.path(output.dir, paste0("Michigan_vs_Eyerich", '_signedpadj_Plot', '.pdf')),
    width = 8, height = 6)
print(
  ggplot2::ggplot(df.Michigan_vs_Eyerich, 
                  aes(x=sign.padj.pso, y=sign.padj.ae, color=diffexpr)) +
    geom_point(size=1) +
    # Add labels
    # Common Upregulated
    geom_label_repel(data = df.Michigan_vs_Eyerich %>% filter(
      sign.padj.pso > 70 & sign.padj.ae < 15 & sign.padj.ae > -log10(0.05)),
      aes(label = hgnc_symbol), max.overlaps = 100, size = 2, color='black',
      nudge_x = 10, nudge_y = -0.15, segment.size = 0.1) +  # label.size = NA -> removes boxes
    # Diff Upregulated down regulated
    geom_label_repel(data = df.Michigan_vs_Eyerich %>% filter(
      sign.padj.ae > 50 & sign.padj.pso < 80 & sign.padj.pso > -log10(0.05)),
      aes(label = hgnc_symbol), max.overlaps = 100, size = 2, color='black',
      nudge_x = -20, nudge_y = 0.15, segment.size = 0.1) +
    geom_label_repel(data = df.Michigan_vs_Eyerich %>% filter(
      sign.padj.ae > 70 & sign.padj.pso > 100),
      aes(label = hgnc_symbol), max.overlaps = 100, size = 2, color='black',
      nudge_x = 20, nudge_y = 1, segment.size = 0.1) +
    # Common Downregulated
    geom_label_repel(data = df.Michigan_vs_Eyerich %>% filter(
      sign.padj.ae < -30 & sign.padj.pso < -10),
      aes(label = hgnc_symbol), max.overlaps = 100, size = 2, color='black',
      nudge_x = -10, nudge_y = -10, segment.size = 0.1) +
    # Choose color
    scale_color_manual(
      name="Differentially regulated",
      labels = c("Down",
                 "Not Sig.",
                 "Up",
                 "Up Eyerich\nDown Michigan",
                 "Up Michigan\nDown Eyerich"),
      values = c('steelblue2','grey', 'indianred2', 'darkblue', 'darkred')) +
    theme_classic() +
    theme(axis.title = element_text(size = axis.title), 
          axis.text = element_text(size = axis.text),
          legend.text = element_text(size = legend.text),
          legend.title = element_text(size = legend.title)) +
    # xlim(-100, 200) +
    # ylim(-50, 100) +
    annotate('text', x = 15, y=-20,
             label = paste0(
               "R = ", format(round(res.Michigan_vs_Eyerich$estimate, 2), nsmall = 2),
               '\np-value = ', format(res.Michigan_vs_Eyerich$p.value, scientific = TRUE))) +
    # geom_smooth(method=lm, se=FALSE, aes(group=1), color='black') +
    labs(x="Signed p-adj. value\nMichigan, L vs NL", 
         y = "Signed p-adj. value\nEyerich, L vs NL",
         color="Differentially regulated")
)
dev.off()





meta.data$age <- ceiling(as.numeric(meta.data$age))
meta.data %>%
  ggplot(aes(x=Sex.x, y=age)) +
  geom_boxplot() + theme_classic() +
  theme(axis.title = element_text(size = 20),axis.text = element_text(size = 15)) +
  labs(x="Sex", y = "Age")

meta.data.new <- as.data.frame(colData(dds))
meta.data.new %>%
  ggplot(aes(x=Sex.x, y=age, color=healthysamp_diag_v2)) +
  geom_boxplot() + theme_classic() +
  theme(axis.title = element_text(size = 20),axis.text = element_text(size = 15)) +
  labs(x="Sex", y = "Age", color='Diagnosis')

meta.data.new %>%
  ggplot(aes(x=Sex.x, y=age, color=Pattern)) +
  geom_boxplot() + theme_classic() +
  theme(axis.title = element_text(size = 20),axis.text = element_text(size = 15)) +
  labs(x="Sex", y = "Age", color='Pattern')

  
  
  