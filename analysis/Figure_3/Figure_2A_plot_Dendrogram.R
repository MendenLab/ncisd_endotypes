rm(list = ls(all=TRUE))

library('ggtree')
library('ggplot2')
library(ggnewscale)
library('factoextra')
library('ape')
library("cowplot")
library("dendextend")

# =========================
#       Parameters
# =========================
path.wd <- '/Volumes/CH__data/Projects/Eyerich_AG_projects/BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer/analysis/Molecular_subtypes'

output_dir <- file.path(path.wd, 'output', 'R_final_Dendrogram_gene_selection')
output_dir <- file.path(output_dir, Sys.Date())
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

date_result_files <- '2023-06-01'

# Plot parameters
fontsize = 5
fontsize_strip = 3.88
offset.text = 0.2
offset.val = 0.3

# HC parameters
opt.method = 'cosine'
opt.linkage_method = 'ward.D2'

# =========================
#       Load data
# =========================
df.grouped.counts <- read.csv(file.path(
  path.wd, 'output', 'Dendrogram', date_result_files, 
  'Normed_Dataframe_mean_gene_selection.csv'))
sorting <- order(as.integer(gsub("E", "\\1", df.grouped.counts$cluster)))
df.grouped.counts <- df.grouped.counts[ sorting, ]
cluster.names <- df.grouped.counts$cluster
row.names(df.grouped.counts) <- seq(1, dim(df.grouped.counts)[1])
df.grouped.counts$cluster <- NULL


# Load probability histograms for Pattern and diag
df.pattern <- read.csv(file.path(path.wd, 'Histogram', date_result_files, 
                                 'Probability_histogram_Pattern.csv'))
df.pattern <- df.pattern[ sorting, ]
row.names(df.pattern) <- seq(1, dim(df.pattern)[1])
df.diag <- read.csv(file.path(path.wd, 'Histogram', date_result_files, 
                              'Probability_histogram_diag.csv'))
df.diag <- df.diag[ sorting, ]
row.names(df.diag) <- seq(1, dim(df.diag)[1])


# =========================
# Hierarchical clustering
# =========================

# Calculate cosine distance
Matrix <- as.matrix(df.grouped.counts)
res.dist <- Matrix / sqrt(rowSums(Matrix * Matrix))
res.dist <- res.dist %*% t(res.dist)
res.dist <- as.dist(1 - res.dist)

# Hierarchical clustering
Z <- stats::hclust(d=res.dist, method=opt.linkage_method)
Z$labels <- cluster.names
tree <- ape::as.phylo(Z)


# =========================
#           Colors
# =========================
# Add majority Pattern and Diagnosis to dendrogram
# First we need to order the Pattern and Diagnosis matrix according to the order 
# in the phylogenetic tree:
# https://yulab-smu.top/treedata-book/chapter7.html
color.vec <- c("#006ddb", "#b6dbff", "#004949", "#009292", "#ff6db6", "#490092",
               "#b66dff", "#000000", "#920000", "#E69F00", "#D55E00", "#8B4513",
               "#999999")
color.vec <- as.factor(color.vec[Z$order])


order_diag_pattern <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13)
df.pattern.colors <- data.frame(Pattern = as.factor(order_diag_pattern))  
rownames(df.pattern.colors) <- tree$tip.label
color.pattern <- c(
  '1'='royalblue', '2'='darkred', '3'='royalblue', '4'='darkorange', '5'='darkred', 
  '6'='grey', '7'='royalblue',  '8'='darkred', '9'='royalblue', '10'='royalblue', 
  '11'='grey', '12'='darkviolet', '13'='darkorange')
df.diag.colors <- data.frame(Diagnosis = as.factor(order_diag_pattern))
rownames(df.diag.colors) <- tree$tip.label
color.diag <- c(
  '1'='blue', '2'='maroon', '3'='blue', '4'='orange', '5'='maroon', '6'='maroon', 
  '7'='blue', '8'='maroon', '9'='blue', '10'='dodgerblue', '11'='azure3', 
  '12'='mediumorchid', '13'='goldenrod')

# Reorder Endotypes in diag and Pattern probability dataframe
df.diag.melted <- reshape2::melt(df.diag, id='Endotypes')
df.diag.melted$Endotypes <- factor(df.diag.melted$Endotypes, 
                                   level=Z$labels[Z$order])
df.pattern.melted <- reshape2::melt(df.pattern, id='Endotypes')
df.pattern.melted$Endotypes <- factor(df.pattern.melted$Endotypes, 
                                      level=Z$labels[Z$order])
df.pattern.melted$variable <- as.factor(df.pattern.melted$variable)


# =========================
#      Plots - circular
# =========================
# Plot the best result as circular tree
pdf(file = file.path(
  output_dir, paste0("Circular_", opt.method, "_",
                     opt.linkage_method, '.pdf')), width = 6, height =6)
p.tree <- ggtree::ggtree(tree, layout='fan') + 
  geom_tiplab(size=fontsize) +
  # geom_hilight( node=15, fill='darkblue', alpha=.5) + 
  # geom_hilight(node=17, fill='darkorange', alpha=.5) + 
  # geom_hilight(node=18, fill='darkred', alpha=.5) + 
  # Highlight three outer nodes
  geom_point2(aes(subset=(node==15)), shape=4, size=3, fill='green', stroke=3) + 
  geom_point2(aes(subset=(node==19)), shape=8, size=3, fill='green', stroke=2) + 
  geom_point2(aes(subset=(node==17)), shape=17, size=4, fill='green') + 
  geom_point2(aes(subset=(node==25)), shape=15, size=4, fill='green') + 
  geom_point2(aes(subset=(node==18)), shape=21, size=4, fill='black') + 
  geom_point2(aes(subset=(node==21)), shape=18, size=4, fill='black') +
  geom_point2(aes(subset=(node==20)), shape=25, size=4, fill='black') 
# # Add names and strips aroiund dendrogramm
# geom_strip('9', '3', barsize=1, color='darkblue', fontsize = fontsize_strip,
#            label = "Pattern 3 like", offset.text=offset.text, offset=offset.val, angle = 60) +
# geom_strip('10', '8', barsize=1, color='darkred', fontsize = fontsize_strip,
#            label = "Pattern 2a like", offset.text=offset.text, offset=offset.val, angle = -200) +
# geom_strip('5', '13', barsize=1, color='darkorange', fontsize = fontsize_strip,
#            label="Pattern 1 like", offset.text=offset.text, offset=offset.val, angle = -75) 
# Add heatmap
p.tree <- gheatmap(p.tree, df.pattern.colors, offset=.2, width=.1,
                   colnames_angle=0, colnames_offset_y = -1) +
  ggplot2::scale_fill_manual(values = color.pattern)
p.tree <- p.tree + ggnewscale::new_scale_fill()
p.tree <- gheatmap(p.tree, df.diag.colors, offset=.32, width=.1,
                   colnames_angle=0, colnames_offset_y = -1) + 
  ggplot2::scale_fill_manual(values = color.diag)
p.tree <- p.tree + theme(legend.position="none") 
# How far to open the circle 190 looks good as well
p.tree <- open_tree(p.tree, 0)  # 
# ggtree::rotate(tree_view = p.tree, node=14)

print(p.tree)
dev.off()


# =========================
#      Plots - rectangular
# =========================
p.tree.rect <- ggtree::ggtree(tree, layout='dendrogram') + 
  geom_tiplab(size=fontsize, nudge_x=-0.05, nudge_y=-0.2) +
  # Highlight three outer nodes
  geom_point2(aes(subset=(node==15)), shape=4, size=4, fill='green', stroke=3) + 
  geom_point2(aes(subset=(node==19)), shape=8, size=4, fill='green', stroke=2) + 
  geom_point2(aes(subset=(node==17)), shape=17, size=5, fill='green') + 
  geom_point2(aes(subset=(node==25)), shape=15, size=5, fill='green') + 
  geom_point2(aes(subset=(node==18)), shape=21, size=5, fill='black') + 
  geom_point2(aes(subset=(node==21)), shape=18, size=5, fill='black') +
  geom_point2(aes(subset=(node==20)), shape=25, size=5, fill='black') +
  # Add names and strips around dendrogramm
  geom_strip('E11', 'E13', barsize=1, color='darkblue', fontsize = fontsize_strip,
             label = "Pattern 3 like", offset.text=0.1, offset=0.1, angle = 0, hjust=0.5) +
  geom_strip('E5', 'E10', barsize=1, color='darkred', fontsize = fontsize_strip,
             label = "Pattern 2a like", offset.text=0.1, offset=0.1, angle = 0, hjust=0.5) +
  geom_strip('E1', 'E4', barsize=1, color='darkorange', fontsize = fontsize_strip,
             label="Pattern 1 like", offset.text=0.1, offset=0.1, angle = 0, hjust=0.5)
# Add heatmap
# p.tree.rect <- ggtree::gheatmap(p.tree.rect, df.pattern.colors, offset=.2, width=.1,
#                colnames_angle=0, colnames_offset_y = -1) +
#   ggplot2::scale_fill_manual(values = color.pattern)
# p.tree.rect <- p.tree.rect + ggnewscale::new_scale_fill()
# p.tree.rect <- ggtree::gheatmap(p.tree.rect, df.diag.colors, offset=.32, width=.1,
#                colnames_angle=0, colnames_offset_y = -1) +
#   ggplot2::scale_fill_manual(values = color.diag)
# p.tree.rect <- p.tree.rect + theme(legend.position="none")
# Add barplot
p.barplot.pattern <- ggplot(df.pattern.melted, aes(fill=variable, y=value, x=Endotypes)) + 
  geom_bar(position="stack", stat="identity", width=0.98) + 
  scale_fill_manual(
    values=c('darkorange', 'darkred', 'tomato', 'royalblue', 'mediumvioletred',
             'magenta', 'darkviolet', 'grey')) +
  theme_classic() + theme(legend.position="none") + xlab(NULL) + ylab(NULL)
# Barplot diag
p.barplot.diag <-  ggplot(df.diag.melted, aes(fill=variable, y=value, x=Endotypes)) + 
  geom_bar(position="stack", stat="identity", width=0.98) + 
  scale_fill_manual(
    values=c('orange', 'goldenrod', 'maroon', 'firebrick', 'salmon', 'blue',
             'dodgerblue', 'palevioletred', 'mediumvioletred', 'deeppink', 'magenta',
             'violet', 'darkviolet', 'mediumorchid', 'dimgrey',
             'darkgray', 'azure3', 'grey', 'slategrey',
             'lightgrey', 'whitesmoke')) +
  theme_classic() + theme(legend.position="none") + xlab(NULL) + ylab(NULL)


pdf(file = file.path(output_dir, paste0("Rectangular_", opt.method, "_",
                                        opt.linkage_method, '.pdf')), 
    width = 8, height =8)
p.tree.rect.draw <- ggdraw() +
  draw_plot(p.tree.rect, x = 0.1, y = 0.42, width = 0.9, height = .6) +
  draw_plot(p.barplot.pattern, x = 0.084, y = .22, width = 0.9, height = .2) +
  draw_plot(p.barplot.diag, x = 0.084, y = 0, width = 0.9, height = 0.2) 
  # draw_plot_label(label = c("A", "B", "C"), size = 15,
  #                 x = c(0, 0., 0), y = c(0.99, 0.66, 0.33))

print(p.tree.rect.draw)
dev.off()


pdf(file = file.path(output_dir, paste0("Final_Rectangular_", opt.method, "_",
                                        opt.linkage_method, '.pdf')), 
    width = 7, height = 4)
p.tree.rect.final <- ggtree::ggtree(tree, layout='dendrogram') + 
  geom_tiplab(size=fontsize, nudge_x=-0.05, nudge_y=-0.2) +
  # Highlight three outer nodes
  geom_point2(aes(subset=(node==14)), shape=23, size=4, fill='white', stroke=2) +
  geom_point2(aes(subset=(node==16)), shape=21, size=4, fill='white', stroke=2) + 
  geom_point2(aes(subset=(node==18)), shape=25, size=5, fill='white', stroke=2) 
print(p.tree.rect.final)
dev.off()



# Print to .pdf
pdf(file = file.path(
  output_dir, paste0("Rectangular_", opt.method, "_",
                     opt.linkage_method, '.pdf')), width = 8, height =4)
print(p.tree.rect)
dev.off()


pdf(file = file.path(
  output_dir, paste0("Histogram_pattern_", opt.method, "_",
                     opt.linkage_method, '.pdf')), width = 8, height=4)
print(p.barplot.pattern)
dev.off()


pdf(file = file.path(
  output_dir, paste0("Histogram_diag_", opt.method, "_", 
                     opt.linkage_method, '.pdf')), width = 8, height=4)
print(p.barplot.diag)
dev.off()

