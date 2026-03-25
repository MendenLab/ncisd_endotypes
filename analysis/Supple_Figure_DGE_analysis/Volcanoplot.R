# library(dplyr)
library(ggplot2)
library("gridExtra")
# library(factoextra)


#' @description Creates a Volcano plot labeling significantly expressed genes and with the option to instead label custom genes
#' 
#' @param df.data : Dataframe with columns log2FoldChange, pvalue, padj
#' @param log2fc.cut : Threshold for log2FC
#' @param add.labels : (optional) Names of predefined gene names
def.get.volcanoplot <- function(df.data, log2fc.cut, padj.cut=0.1, 
                                add.labels=NULL, add.annotations=TRUE) 
{
  # 1.Mark diff expressed genes
  df.data <- df.data[order(df.data$padj), ]
  # add a column of NAs
  df.data$diffexpressed <- "Not Sig"
  # if log2Foldchange > 1 and padj < 0.05, set as "Up" 
  df.data$diffexpressed[
    df.data$log2FoldChange > log2fc.cut & df.data$padj < padj.cut] <- "Up"
  # if log2Foldchange < -1 and padj < 0.05, set as "Down"
  df.data$diffexpressed[
    df.data$log2FoldChange < -log2fc.cut & df.data$padj < padj.cut] <- "Down"
  df.data$diffexpressed <- as.factor(df.data$diffexpressed)
  df.data$diffexpressed <- factor(
    df.data$diffexpressed , levels = c("Not Sig", "Up", "Down"))
  
  # 2. Add labels to genes
  df.data <- def.getlabels(func_df.data=df.data, add.labels=NULL) 
  
  # 2.1 DOWN: split into two columns if label length > 20
  mask.down <- !is.na(df.data$label) & df.data$diffexpressed == "Down"
  # mask.down <- mask.down | df.data$label %in% add.labels
  ind.temp.down <- rownames(df.data[mask.down, ])
  second.col.down <- df.data$label[mask.down]
  second.half.down <- ceiling(length(second.col.down) / 2)
  if (length(second.col.down) %% 2 == 0 & length(second.col.down) > 1) 
  {
    second.half.down <- second.half.down + 1
  }
  first.half.down <- length(second.col.down) - second.half.down
  
  # 2.2 UP: split into two columns if label length > 20
  mask.up <- !is.na(df.data$label) & df.data$diffexpressed == "Up"
  # mask.up <- mask.up | df.data$label %in% add.labels 
  ind.temp.up <-  rownames(df.data[mask.up, ])
  second.col.up <- df.data$label[mask.up]
  second.half.up <- ceiling(length(second.col.up) / 2)
  first.half.up <- length(second.col.up) - second.half.up
  if (length(second.col.up) %% 2 == 0 & length(second.col.up) > 1) 
  {
    second.half.up <- second.half.up + 1
  }
  
  
  # 2.3 responder genes - exclude them from labels and set to NaN again
  df.data$resp_label <- NA
  # Add custom labels -> Todo check in CEBPB
  ind.labels <- match(add.labels, df.data$hgnc_symbol)
  # remove NAs
  ind.labels <- ind.labels[!is.na(ind.labels)]
  # in label set again to NA
  df.data$label[ind.labels] <- NA
  # add labels to resp_label
  df.data$resp_label[ind.labels] <- df.data$hgnc_symbol[ind.labels]
  
  ##   split by up and down:
  # 2.3.1 Down responder genes
  mask.resps.down <- df.data$resp_label %in% add.labels & df.data$diffexpressed == "Down"
  ind.temp.resp.down <-  rownames(df.data[mask.resps.down, ])
  second.col.resps.down <- df.data$resp_label[mask.resps.down]
  second.half.resps.down <- ceiling(length(second.col.resps.down) / 2)
  first.half.resps.down <- length(second.col.resps.down) - second.half.resps.down
  # To avoid that labels occur twice if second.col.resps.down is even
  if (length(second.col.resps.down) %% 2 == 0 & length(second.col.resps.down) > 1) 
  {
    second.half.resps.down <- second.half.resps.down + 1
  }
  
  # 2.3.2 UP responder genes
  mask.resps.up <- df.data$resp_label %in% add.labels & df.data$diffexpressed == "Up"
  ind.temp.resp.up <-  rownames(df.data[mask.resps.up, ])
  second.col.resps.up <- df.data$resp_label[mask.resps.up]
  second.half.resps.up <- ceiling(length(second.col.resps.up) / 2)
  first.half.resps.up <- length(second.col.resps.up) - second.half.resps.up
  # To avoid that labels occur twice if second.col.resps.down is even
  if (length(second.col.resps.up) %% 2 == 0 & length(second.col.resps.up) > 1) 
  {
    second.half.resps.up <- second.half.resps.up + 1
  }
  
  # 2.3.3 inbetween responder genes, negative log2FC
  mask.resps.inbe.neg <- df.data$resp_label %in% add.labels & 
    df.data$diffexpressed == "Not Sig" & df.data$log2FoldChange < 0
  ind.temp.resp.inbe.neg <-  rownames(df.data[mask.resps.inbe.neg, ])
  
  # 2.3.3 inbetween responder genes, negative log2FC
  mask.resps.inbe.pos <- df.data$resp_label %in% add.labels & 
    df.data$diffexpressed == "Not Sig" & df.data$log2FoldChange > 0
  ind.temp.resp.inbe.pos <-  rownames(df.data[mask.resps.inbe.pos, ])
  
  
  # 3. set parameters
  # 3.1 find xaxis max limit
  xmax = max(abs(df.data$log2FoldChange[
    !is.na(df.data$log2FoldChange) & (!is.infinite(df.data$log2FoldChange))]))
  xlim_magnitude = floor(log10(xmax))
  # 3.2 Get second largest y-value
  y = -log10(df.data$padj)
  ymax = max(y)
  ymax_second = sort(y)[length(y) - 1]  
  ylim_magnitude = floor(log10(max(-log10(df.data$padj))))
  # 3.3 Point and label size
  size.annot <- 4  # label fontsize
  point.size <- 1
  fill <- NA  # no background color of label boxex
  fontsize <- 1
  segment.size <- 0.3   # size of line pointing from label to point
  box.padding <- 0.1
  color_up <- 'black' # 'darkred'
  color_down <- 'black' # 'darkblue'
  
  # 4. Create Volcano plot
  func.volcanoplot <- ggplot2::ggplot(
    data=df.data %>% plyr::arrange(diffexpressed), 
    ggplot2::aes(x=log2FoldChange, y=-log10(padj), color=diffexpressed, label=label)) +
    ggplot2::geom_point(size=point.size) + 
    ggplot2::theme_classic()
  
  if (add.annotations == TRUE) 
  {
    func.volcanoplot <- func.volcanoplot + ggplot2::scale_x_continuous(
      limits = c(-xmax - 12.5 * 10**xlim_magnitude, xmax + 12.5 * 10**xlim_magnitude))
  } else 
  {
    func.volcanoplot <- func.volcanoplot + ggplot2::scale_x_continuous(
      limits = c(-xmax - 1 * 10**xlim_magnitude, xmax + 1 * 10**xlim_magnitude)) + 
      ggplot2::scale_y_continuous(
        limits = c(0, ymax + 3 * 10**xlim_magnitude))
  }
  
  # 4.1 Split data into two columns
  # Annotate resp genes
  # Down regulated responder genes
  if (length(ind.temp.resp.down) > 0) 
  {
    func.volcanoplot <- func.volcanoplot + ggrepel::geom_label_repel(
      data=df.data[ind.temp.resp.down[1:first.half.resps.down], ], 
      ggplot2::aes(label=df.data[ind.temp.resp.down[1:first.half.resps.down], 'resp_label'],
                   fontsize=fontsize, size=6),
      fill = fill, color = 'darkgoldenrod',
      xlim = c(-xmax - 6*10**xlim_magnitude, -xmax - 6.4*10**xlim_magnitude),
      ylim = c(min(-log10(df.data$padj)), max(-log10(df.data$padj))),
      size=size.annot,  max.overlaps = Inf, 
      direction    = "y",
      hjust        = 0,
      segment.size = segment.size,
      label.size = NA,  # remove boxes around labels
      size = fontsize, 
      box.padding  = box.padding, 
      show.legend  = FALSE) 
    if (length(second.col.resps.down) > 1) 
    {
      func.volcanoplot <- func.volcanoplot + ggrepel::geom_label_repel(
        data=df.data[ind.temp.resp.down[second.half.resps.down:length(second.col.resps.down)], ], 
        ggplot2::aes(label=df.data[
          ind.temp.resp.down[second.half.resps.down:length(second.col.resps.down)], 'resp_label'],
          fontsize=fontsize, size=fontsize), 
        fill = fill, color = 'darkgoldenrod',
        xlim = c(-xmax - 8.4*10**xlim_magnitude, -xmax - 8.8*10**xlim_magnitude),
        ylim = c(min(-log10(df.data$padj)), max(-log10(df.data$padj))),
        size=size.annot,  max.overlaps = Inf, 
        direction    = "y",
        hjust        = 0,
        segment.size = segment.size,
        label.size = NA,  # remove boxes around labels
        box.padding  = box.padding, 
        show.legend  = FALSE)
    } 
  }
  if (length(ind.temp.resp.up) > 0) 
  {
    # UP regulated responder genes
    func.volcanoplot <- func.volcanoplot + ggrepel::geom_label_repel(
      data=df.data[ind.temp.resp.up[1:first.half.resps.up], ],
      ggplot2::aes(label=df.data[ind.temp.resp.up[1:first.half.resps.up], 'resp_label'], 
                   size=fontsize),
      fill = fill, color = 'darkgoldenrod',
      xlim = c(xmax + 6*10**xlim_magnitude - 4*10**(xlim_magnitude-1), xmax + 6.4*10**xlim_magnitude),
      ylim = c(min(-log10(df.data$padj)), max(-log10(df.data$padj))),
      size=size.annot,  max.overlaps = Inf, 
      direction    = "y",
      hjust        = 0,
      segment.size = segment.size,
      label.size = NA,  # remove boxes around labels
      size = fontsize,
      box.padding  = box.padding, show.legend  = FALSE)
    if (length(second.col.resps.up) > 1) 
    {
      func.volcanoplot <- func.volcanoplot + ggrepel::geom_label_repel(
        data=df.data[ind.temp.resp.up[second.half.resps.up:length(second.col.resps.up)], ], 
        ggplot2::aes(label=df.data[ind.temp.resp.up[second.half.resps.up:length(second.col.resps.up)],
                                   'resp_label'], fontsize=5), 
        fill = fill, color = 'darkgoldenrod',
        xlim = c(xmax + 8.4*10**xlim_magnitude - 4*10**(xlim_magnitude-1), 
                 xmax + 8.8*10**xlim_magnitude),
        ylim = c(min(-log10(df.data$padj)), max(-log10(df.data$padj))),
        size=size.annot,  max.overlaps = Inf, direction    = "y", size = fontsize, 
        label.size = NA,  # remove boxes around labels
        hjust        = 0, segment.size = segment.size, box.padding  = box.padding, cex.lab=1.5,
        show.legend  = FALSE) 
    } 
  }
  if (length(ind.temp.resp.inbe.neg) > 0) 
  {
    # Inbetween Not sig responder genes
    func.volcanoplot <- func.volcanoplot + ggrepel::geom_label_repel(
      data=df.data[ind.temp.resp.inbe.neg, ],
      ggplot2::aes(label=df.data[ind.temp.resp.inbe.neg, 'resp_label'], 
                   fontsize=fontsize, size=fontsize),
      fill = fill, color = 'grey',
      xlim = c(-xmax,  -xmax + 0.01*10**xlim_magnitude),
      # ylim = c(ymax_second + 3*10**(ylim_magnitude -1), ymax - 3*10**(ylim_magnitude -1)),
      size=size.annot,  max.overlaps = Inf, 
      direction    = "y",
      hjust        = 0,
      segment.size = segment.size,
      size = fontsize, 
      box.padding  = box.padding, cex.lab=1.5, show.legend  = FALSE) 
  }
  if (length(ind.temp.resp.inbe.pos) > 0) 
  {
    func.volcanoplot <- func.volcanoplot + ggrepel::geom_label_repel(
      data=df.data[ind.temp.resp.inbe.pos, ], 
      ggplot2::aes(label=df.data[ind.temp.resp.inbe.pos, 'resp_label'], 
                   fontsize=fontsize, size=fontsize), 
      fill = fill, color = 'grey',
      xlim = c(xmax, xmax + 0.01*10**xlim_magnitude),
      # ylim = c(ymax_second + 3*10**(ylim_magnitude -1), ymax - 3*10**(ylim_magnitude -1)),
      size=size.annot,  max.overlaps = Inf, direction    = "y", size = fontsize, 
      hjust        = 0, segment.size = segment.size, box.padding  = box.padding, cex.lab=1.5,
      show.legend  = FALSE)
  }
  
  
  # Split into 20 each; use different nudges
  # Anotate top diff expressed genes
  # Donw regulated genes
  if (add.annotations == TRUE) 
  {
    func.volcanoplot <- func.volcanoplot + ggrepel::geom_label_repel(
      data=df.data[ind.temp.down[1:first.half.down], ], 
      fill = fill,  color = color_down,
      xlim = c(-xmax - 14.*10**xlim_magnitude + 14*10**(xlim_magnitude-1), 
               -xmax - 14.4*10**xlim_magnitude),
      # xlim = c(-xmax - 4.*10**xlim_magnitude + 4*10**(xlim_magnitude-1), 
      #          -xmax - 4.4*10**xlim_magnitude),
      # xlim = c(NA, -xmax - 4.4*10**xlim_magnitude), 
      ylim = c(min(-log10(df.data$padj)), max(-log10(df.data$padj))),
      size=size.annot,  max.overlaps = Inf, 
      direction    = "y",
      hjust        = 0,
      segment.size = segment.size,
      label.size = NA,  # remove boxes around labels
      size = fontsize, aes(fontsize=fontsize, size=fontsize),
      box.padding  = box.padding, cex.lab=1.5, show.legend  = FALSE)
    if (length(second.col.down) > 1) 
    {
      func.volcanoplot <- func.volcanoplot + ggrepel::geom_label_repel(
        data=df.data[ind.temp.down[second.half.down:length(second.col.down)], ], 
        fill = fill, color = color_down,
        xlim = c(-xmax - 6.4*10**xlim_magnitude + 10*10**(xlim_magnitude-1),
                 -xmax - 6.1*10**xlim_magnitude),
        # xlim = c(-xmax - 0.4*10**xlim_magnitude + 4*10**(xlim_magnitude-1),
        #          -xmax - 0.1*10**xlim_magnitude),
        # xlim = c(NA, -xmax - 0.1*10**xlim_magnitude), 
        ylim = c(min(-log10(df.data$padj)), max(-log10(df.data$padj))),
        size=size.annot,  max.overlaps = Inf, 
        direction    = "y",
        hjust        = 0,
        segment.size = segment.size,
        label.size = NA,  # remove boxes around labels
        size = fontsize, aes(fontsize=fontsize, size=fontsize), 
        box.padding  = box.padding, cex.lab=1.5, show.legend  = FALSE)
    } 
    
    # Up regulated genes
    func.volcanoplot <- func.volcanoplot + ggrepel::geom_label_repel(
      data=df.data[ind.temp.up[1:first.half.up], ], 
      fill = fill,  color = color_up,
      xlim = c(xmax + 6.1*10**xlim_magnitude - 10*10**(xlim_magnitude-1),
               xmax + 6.4*10**xlim_magnitude),
      # xlim = c(xmax + 0.1*10**xlim_magnitude - 4*10**(xlim_magnitude-1),
      #          xmax + 0.4*10**xlim_magnitude),
      # xlim = c(xmax + 0.1*10**xlim_magnitude - 4*10**(xlim_magnitude-1), NA), 
      size=size.annot,  max.overlaps = Inf,
      direction    = "y",
      hjust        = 0,
      segment.size = segment.size,
      label.size = NA,  # remove boxes around labels
      size = fontsize, aes(fontsize=fontsize, size=fontsize),
      box.padding  = box.padding, cex.lab=1.5, show.legend  = FALSE)
    if (length(second.col.up) > 1) 
    {
      func.volcanoplot <- func.volcanoplot +  ggrepel::geom_label_repel(
        data=df.data[ind.temp.up[second.half.up:length(second.col.up)], ], 
        fill = fill,  color = color_up,
        xlim = c(xmax + 14*10**xlim_magnitude - 14*10**(xlim_magnitude-1),
                 xmax + 14.4*10**xlim_magnitude),
        # xlim = c(xmax + 4*10**xlim_magnitude - 4*10**(xlim_magnitude-1),
        #          xmax + 4.4*10**xlim_magnitude),
        # xlim = c(xmax + 4*10**xlim_magnitude - 4*10**(xlim_magnitude-1), NA), 
        size=size.annot,  max.overlaps = Inf,
        direction    = "y",
        hjust        = 0,
        segment.size = segment.size,
        label.size = NA,  # remove boxes around labels
        size = fontsize, aes(fontsize=fontsize, size=fontsize),
        box.padding  = box.padding, show.legend  = FALSE)
    }
    
  }
  
  # Legend 
  if ((length(ind.temp.up) >= 1) & (length(ind.temp.down) < 1)) 
  {
    func.volcanoplot <- func.volcanoplot + 
      ggplot2::labs(color = "Diff. regulated", 
                    x = bquote(log[2] ~ "FC"), y=bquote(log[10] ~ "p-adj value")) +
      ggplot2::scale_color_manual(
        values=c('Not Sig' = "grey", "Up" = "indianred2"), 
        labels = c("Not Sig.", "Up")) +
      ggplot2::theme(axis.title = element_text(size = 18), axis.text = element_text(size = 14),
                     legend.text = element_text(size = 14), legend.title = element_text(size = 16))  
  } else if ((length(ind.temp.up) < 1) & (length(ind.temp.down) >= 1)) 
  {
    func.volcanoplot <- func.volcanoplot + 
      ggplot2::labs(color = "Diff. regulated", 
                    x = bquote(log[2] ~ "FC"), y=bquote(log[10] ~ "p-adj value")) +
      ggplot2::scale_color_manual(
        values=c('Not Sig' = "grey", 'Down' = "steelblue2"), 
        labels = c("Not Sig.", "Down")) +
      ggplot2::theme(axis.title = element_text(size = 18), axis.text = element_text(size = 14),
                     legend.text = element_text(size = 14), legend.title = element_text(size = 16)) 
  } else 
  {
    func.volcanoplot <- func.volcanoplot + 
      ggplot2::labs(color = "Diff. regulated", 
                    x = bquote(log[2] ~ "FC"), y=bquote(log[10] ~ "p-adj value")) +
      ggplot2::scale_color_manual(
        values=c('Not Sig' = "grey", 'Down' = "steelblue2", "Up" = "indianred2"), 
        labels = c("Not Sig.", "Up", "Down")) +
      ggplot2::theme(text = element_text(family = "sans"),
                     axis.title = element_text(size = 18), 
                     axis.text = element_text(size = 14),
                     legend.text = element_text(size = 14), 
                     legend.title = element_text(size = 16))  
  }

  
  return(func.volcanoplot)
}


def.getlabels <- function(func_df.data, add.labels) 
{
  if (!is.null(add.labels)) 
  {
    func_df.data$label <- NA
    # Add custom labels -> Todo check in CEBPB
    ind.labels <- match(add.labels, func_df.data$hgnc_symbol)
    # remove NAs
    ind.labels <- ind.labels[!is.na(ind.labels)]
    func_df.data$label[ind.labels] <- func_df.data$hgnc_symbol[ind.labels]
    
    # add top 10 DEGs annotations
    top_genes <- 10
    mask.sig <- func_df.data$diffexpressed == "Up"
    max_ele <- top_genes + 1  # sum(mask.sig, na.rm = TRUE)
    func_df.data$label[mask.sig][1:max_ele] <- func_df.data$hgnc_symbol[mask.sig][1:max_ele]
    
    mask.sig <- func_df.data$diffexpressed == "Down"
    func_df.data$label[mask.sig][1:max_ele] <- func_df.data$hgnc_symbol[mask.sig][1:max_ele]
    
  } else 
  {
    # Create a new column "label", contains the name of DEx genes (NA in case they are not)
    func_df.data$label <- NA
    
    # add top 10 DEGs annotations
    top_genes <- 10
    mask.sig <- func_df.data$diffexpressed == "Up"
    max_ele <- top_genes + 1  # sum(mask.sig, na.rm = TRUE)
    func_df.data$label[mask.sig][1:max_ele] <- func_df.data$hgnc_symbol[mask.sig][1:max_ele]
    
    mask.sig <- func_df.data$diffexpressed == "Down"
    func_df.data$label[mask.sig][1:max_ele] <- func_df.data$hgnc_symbol[mask.sig][1:max_ele]
  }
  
  # sort by label 
  func_df.data <- func_df.data[order(func_df.data$label), ]
  return(func_df.data)
}
