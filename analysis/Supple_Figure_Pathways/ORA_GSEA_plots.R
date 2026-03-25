# library(forcats)
# library(ggplot2)

def.plot.barplot <- function(enrich_gsea.result, method, showCategories, title, 
                             width, height, output.dir, func_colours=NULL) 
{
  if (typeof(showCategories) != 'character') 
  {
    if (showCategories > nrow(enrich_gsea.result)) 
    {
      showCategories = nrow(enrich_gsea.result)
    }
  } else 
  {
    # get intersection Pathways which are in enrichObject and in the user provided Pathway list
    showCategories = intersect(enrich_gsea.result$Description, showCategories)
  }
  
  pdf(file = file.path(output.dir, title), width = width, height = height)
  # p1 <- ggplot(enrich_gsea.result, showCategory=showCategories, 
  #              aes(NES, forcats::fct_reorder(Description, NES), fill=p.adjust)) + 
  #   geom_col() + theme_classic() + 
  #   ggtitle(method) +
  #   xlab("Normalized Enrichment Score (NES)") +
  #   ylab(NULL) 
  # 
  # if (!is.null(func_colours))
  # {
  #   p1 <- p1 + scale_fill_gradientn(colours=c("#b3eebe", "#46bac2", "#371ea3"),
  #                                   guide=guide_colorbar(reverse=TRUE)) 
  # }
  # 
  
  max.val <- ceiling(max(enrich_gsea.result@result$NES)) + 3
  min.val <- floor(min(enrich_gsea.result@result$NES)) - 3
  
  p1 <- ggplot2::ggplot(
    enrich_gsea.result, showCategory=showCategories,
    ggplot2::aes(x = NES, y = forcats::fct_reorder(Description, NES), fill=p.adjust)) +
    ggplot2::geom_bar(stat = "identity",
             show.legend = TRUE,
             ggplot2::aes(fill = p.adjust),     
             color = "white") +
    ggplot2::geom_hline(yintercept = 0, color = 1, lwd = 0.2) +
    ggplot2::geom_text(
      ggplot2::aes(label = Description, # Text with groups
                   hjust = ifelse(NES < 0, 1., 0),
                   vjust = 0.5), size = 2) +
    ggplot2::xlab("Normalized Enrichment Score (NES)") +
    ggplot2::ylab(NULL) +
    ggplot2::scale_x_continuous(breaks = seq(min.val, max.val, by = 1),
                       limits = c(min.val, max.val)) +
    ggplot2::scale_fill_gradient2(low = "lightblue",
                         mid = "aliceblue",
                         high = "pink") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),  # Remove Y-axis texts
          axis.ticks.y = ggplot2::element_blank(), # Remove Y-axis ticks
          panel.grid.major.y = ggplot2::element_blank(), # Remove horizontal grid
          panel.grid.major.x = ggplot2::element_blank(), # Remove vertical grid
          # panel.grid.minor = element_blank(),
          # panel.grid.major = element_blank(), 
          panel.background = ggplot2::element_blank() )

  
  cowplot::plot_grid(p1, ncol = 1, nrow = 1)
  print(p1)
  dev.off()
}


def.group.dotplot <- function(enrich_gsea.result, method, showCategories, title, 
                              width, height, output.dir, groups, func_colours=NULL) 
{
  if (typeof(showCategories) != 'character') 
  {
    if (showCategories > nrow(enrich_gsea.result)) 
    {
      showCategories = nrow(enrich_gsea.result)
    }
  } else 
  {
    # get intersection Pathways which are in enrichObject and in the user provided Pathway list
    showCategories = intersect(enrich_gsea.result$Description, showCategories)
  }
  
  pdf(file = file.path(output.dir, title), width = width, height = height)
  p1 <- enrichplot::dotplot(
    enrich_gsea.result, showCategory = showCategories, 
    x = "GeneRatio", color = "p.adjust", 
    title = paste("Enriched ", method , "\nPathways of ", groups, sep="") , 
    split=".sign") + 
    facet_grid(.~.sign) + 
    ggplot2::xlab("GeneRatio") +
    if (!is.null(func_colours))
    {
      p1 <- p1 + scale_fill_gradientn(colours=c("#b3eebe", "#46bac2", "#371ea3"),
                                      guide=guide_colorbar(reverse=TRUE)) 
    }
  
  cowplot::plot_grid(p1, ncol = 1, nrow = 1)
  print(p1)
  dev.off()
}



def.group.barplot <- function(enrich_gsea.result, method, showCategories, title, 
                              width, height, output.dir, groups, func_colours=NULL) 
{
  if (typeof(showCategories) != 'character') 
  {
    if (showCategories > nrow(enrich_gsea.result)) 
    {
      showCategories = nrow(enrich_gsea.result)
    }
  } else 
  {
    # get intersection Pathways which are in enrichObject and in the user provided Pathway list
    showCategories = intersect(enrich_gsea.result$Description, showCategories)
  }
  
  pdf(file = file.path(output.dir, title), width = width, height = height)
  p1 <- enrichplot::dotplot(enrich_gsea.result, showCategory = showCategories, 
                            title = paste("Enriched ", method , " Pathways of ", groups, sep="") , 
                            split=".sign") + 
    facet_grid(.~.sign) + 
    xlab("GeneRatio") +
    if (!is.null(func_colours))
    {
      p1 <- p1 + scale_fill_gradientn(colours=c("#b3eebe", "#46bac2", "#371ea3"),
                                      guide=guide_colorbar(reverse=TRUE)) 
    }
  
  cowplot::plot_grid(p1, ncol = 1, nrow = 1)
  print(p1)
  dev.off()
}

def.pathways.cnetplot <- function(enrich.result, method, entrezid_log2fc, 
                                  showCategories, title, width, height, output.dir) 
{
  if (typeof(showCategories) != 'character') 
  {
    if (showCategories > nrow(enrich.result)) 
    {
      showCategories = nrow(enrich.result)
    }
  } else 
  {
    # get intersection Pathways which are in enrichObject and in the user provided Pathway list
    showCategories = intersect(enrich.result$Description, showCategories)
  }
  
  pdf(file = file.path(output.dir, title), width = width, height = height)
  p1 = enrichplot::cnetplot(enrich.result, showCategory = showCategories, 
                            categorySize = "p.adjust", color.params = list(
                              edge = TRUE, foldChange = entrezid_log2fc), max.overlaps=Inf) + 
    ggplot2::ggtitle(paste(method, sep = " ")) 
  
  cowplot::plot_grid(p1, ncol = 1, nrow = 1)
  print(p1)
  dev.off()
}



def.plot.enricher_barplot <- function(enricher.result, method, showCategories, title, 
                             width, height, output.dir, func_colours=NULL) 
{
  if (typeof(showCategories) != 'character') 
  {
    if (showCategories > nrow(enricher.result)) 
    {
      showCategories = nrow(enricher.result)
    }
  } else 
  {
    # get intersection Pathways which are in enrichObject and in the user provided Pathway list
    showCategories = intersect(enricher.result$Description, showCategories)
  }
  
  pdf(file = file.path(output.dir, title), width = width, height = height)
  p1 = barplot(enricher.result, showCategory = showCategories, 
               color='p.adjust', x='GeneRatio')  +
    ggplot2::xlab("GeneRatio") + ggplot2::ylab('Pathways') + 
    ggplot2::theme_classic() + ggplot2::ggtitle(method)
  if (!is.null(func_colours))
  {
    p1 <- p1 + scale_fill_gradientn(colours=c("#b3eebe", "#46bac2", "#371ea3"),
                                    guide=guide_colorbar(reverse=TRUE)) 
  }
  
  cowplot::plot_grid(p1, ncol = 1, nrow = 1)
  print(p1)
  dev.off()
}
