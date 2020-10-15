plot_spatial_GO_term <- function(pathway_add, plot_tumor = T, plot_interface = F) {
  # this assumes that the Seurat objects for ABC are already loaded, and the expression and pval matrices
  
  if (length(pathway.add) == 1) {
    
    pathway.add <- pathway_add
    samples <- c("A1", "B1", "C1")
    plots <- NULL
    
    if (plot_tumor) {
      for (sample in samples) {
        if (sample == "A1") {
          metadata.add <- A1.t.expr[rownames(A1.t.expr) %in% pathway.add,] %>% t() %>% rescale(to = c(-1,1)) %>% data.frame() %>% rownames_to_column() %>% deframe()
          names(metadata.add) <- gsub(pattern = "_1", replacement = "", x = names(metadata.add))
          A1 <- AddMetaData(A1, metadata = metadata.add, col.name = pathway.add)
          pval <- A1.t.pval[rownames(A1.t.pval) == pathway.add,] %>% signif(digits = 4)
          plots[[sample]] <- SpatialPlot(A1, features = pathway.add, image.alpha = 0, stroke = 0, pt.size.factor = 1.3) + 
            scale_fill_gradientn(colours = brewer.rdbu(n = 100) %>% rev()) +
            labs(title = pathway.add, subtitle = paste("P =", pval)) +
            theme(legend.position = "none", plot.title = element_text(size = 14, hjust = 0.5), plot.subtitle = element_text(size = 12, hjust = 0.5))
        } else if (sample == "B1") {
          metadata.add <- B1.t.expr[rownames(B1.t.expr) %in% pathway.add,] %>% t() %>% rescale(to = c(-1,1)) %>% data.frame() %>% rownames_to_column() %>% deframe()
          names(metadata.add) <- gsub(pattern = "_2", replacement = "", x = names(metadata.add))
          B1 <- AddMetaData(B1, metadata = metadata.add, col.name = pathway.add)
          pval <- B1.t.pval[rownames(B1.t.pval) == pathway.add,] %>% signif(digits = 4)
          plots[[sample]] <- SpatialPlot(B1, features = pathway.add, image.alpha = 0, stroke = 0, pt.size.factor = 1.3) + 
            scale_fill_gradientn(colours = brewer.rdbu(n = 100) %>% rev()) +
            labs(title = pathway.add, subtitle = paste("P =", pval)) +
            theme(legend.position = "none", plot.title = element_text(size = 14, hjust = 0.5), plot.subtitle = element_text(size = 12, hjust = 0.5))
        } else if (sample == "C1") {
          metadata.add <- C1.t.expr[rownames(C1.t.expr) %in% pathway.add,] %>% t() %>% rescale(to = c(-1,1)) %>% data.frame() %>% rownames_to_column() %>% deframe()
          names(metadata.add) <- gsub(pattern = "_3", replacement = "", x = names(metadata.add))
          C1 <- AddMetaData(C1, metadata = metadata.add, col.name = pathway.add)
          pval <- C1.t.pval[rownames(C1.t.pval) == pathway.add,] %>% signif(digits = 4)
          plots[[sample]] <- SpatialPlot(C1, features = pathway.add, image.alpha = 0, stroke = 0, pt.size.factor = 1.3) + 
            scale_fill_gradientn(colours = brewer.rdbu(n = 100) %>% rev()) +
            labs(title = pathway.add, subtitle = paste("P =", pval)) +
            theme(legend.position = "none", plot.title = element_text(size = 14, hjust = 0.5), plot.subtitle = element_text(size = 12, hjust = 0.5))
        }
      }
    } else if (plot_interface) {
      for (sample in samples) {
        if (sample == "A1") {
          metadata.add <- A1.i.expr[rownames(A1.i.expr) %in% pathway.add,] %>% t() %>% rescale(to = c(-1,1)) %>% data.frame() %>% rownames_to_column() %>% deframe()
          names(metadata.add) <- gsub(pattern = "_1", replacement = "", x = names(metadata.add))
          A1 <- AddMetaData(A1, metadata = metadata.add, col.name = pathway.add)
          pval <- A1.i.pval[rownames(A1.i.pval) == pathway.add,] %>% signif(digits = 4)
          plots[[sample]] <- SpatialPlot(A1, features = pathway.add, image.alpha = 0, stroke = 0, pt.size.factor = 1.3) + 
            scale_fill_gradientn(colours = brewer.rdbu(n = 100) %>% rev()) +
            labs(title = pathway.add, subtitle = paste("P =", pval)) +
            theme(legend.position = "none", plot.title = element_text(size = 14, hjust = 0.5), plot.subtitle = element_text(size = 12, hjust = 0.5))
        } else if (sample == "B1") {
          metadata.add <- B1.i.expr[rownames(B1.i.expr) %in% pathway.add,] %>% t() %>% rescale(to = c(-1,1)) %>% data.frame() %>% rownames_to_column() %>% deframe()
          names(metadata.add) <- gsub(pattern = "_2", replacement = "", x = names(metadata.add))
          B1 <- AddMetaData(B1, metadata = metadata.add, col.name = pathway.add)
          pval <- B1.i.pval[rownames(B1.i.pval) == pathway.add,] %>% signif(digits = 4)
          plots[[sample]] <- SpatialPlot(B1, features = pathway.add, image.alpha = 0, stroke = 0, pt.size.factor = 1.3) + 
            scale_fill_gradientn(colours = brewer.rdbu(n = 100) %>% rev()) +
            labs(title = pathway.add, subtitle = paste("P =", pval)) +
            theme(legend.position = "none", plot.title = element_text(size = 14, hjust = 0.5), plot.subtitle = element_text(size = 12, hjust = 0.5))
        } else if (sample == "C1") {
          metadata.add <- C1.i.expr[rownames(C1.i.expr) %in% pathway.add,] %>% t() %>% rescale(to = c(-1,1)) %>% data.frame() %>% rownames_to_column() %>% deframe()
          names(metadata.add) <- gsub(pattern = "_3", replacement = "", x = names(metadata.add))
          C1 <- AddMetaData(C1, metadata = metadata.add, col.name = pathway.add)
          pval <- C1.i.pval[rownames(C1.i.pval) == pathway.add,] %>% signif(digits = 4)
          plots[[sample]] <- SpatialPlot(C1, features = pathway.add, image.alpha = 0, stroke = 0, pt.size.factor = 1.3) + 
            scale_fill_gradientn(colours = brewer.rdbu(n = 100) %>% rev()) +
            labs(title = pathway.add, subtitle = paste("P =", pval)) +
            theme(legend.position = "none", plot.title = element_text(size = 14, hjust = 0.5), plot.subtitle = element_text(size = 12, hjust = 0.5))
        }
      }
    }
    spatial_plots <- plot_grid(plotlist = plots, nrow = 1)
    return(spatial_plots)
  }
}

