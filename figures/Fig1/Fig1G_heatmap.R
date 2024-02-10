library(Seurat)
library(glue)
library(Nebulosa)
library(tidyverse)

scrna_to_plot <- function(scrna){
  if(ncol(scrna)<=20000){
    return(scrna)
  }
  return(subset(scrna, cells=sample(colnames(scrna))[1:20000]))
}



source("../util.R")


object <- load_object("../save/scrna_phase_comparing.Rds")



object$annotation <- factor(object$annotation, levels=c(glue("PT{c(1,2,3,5,6,9)}"), "eGFP+", "STCs"))




# DE plots
cluster_de <- object@tools[["de_annotation"]]
cluster_de <- cluster_de[sapply(cluster_de, function(m) nrow(m) > 0)]

cluster_de_top10 <- lapply(cluster_de, function(x) {
    if("avg_logFC" %in% names(x)){ ## compatible for seurat3
      x$avg_log2FC <- x$avg_logFC / log(2)
    }
    x %>% top_n(10, avg_log2FC) %>% arrange(-avg_log2FC)
})

#cluster_de_top10 <- cluster_de_top10

## top10 DE heatmaps
genes <- as.vector(unlist(sapply(cluster_de_top10, function(x)x$gene)))
object <- ScaleData(object, rownames(object))
col_def <- rev(ggsci::pal_simpsons()(length(unique(object@meta.data[, "annotation"]))))



pdf("save/fig1G.pdf")
plthm = DoHeatmap(
    ## >30k will failed to plot, here we subset when cell number > 20,000
    object = scrna_to_plot(object),
    features = genes,
    group.by = "annotation",
    group.colors = col_def,
    disp.min = -2,
    disp.max = 2,
    slot = "scale.data",
    assay = "RNA",
    raster = FALSE,
    combine = TRUE
  )
print(plthm)
dev.off()


pdf("save/fig1G_raster.pdf")
plthm = DoHeatmap(
    ## >30k will failed to plot, here we subset when cell number > 20,000
    object = scrna_to_plot(object),
    features = genes,
    group.by = "annotation",
    group.colors = col_def,
    disp.min = -2,
    disp.max = 2,
    slot = "scale.data",
    assay = "RNA",
    raster = TRUE,
    combine = TRUE
  )
print(plthm)

plthm = DoHeatmap(
    ## >30k will failed to plot, here we subset when cell number > 20,000
    object = scrna_to_plot(object),
    features = genes,
    group.by = "annotation",
    group.colors = col_def,
    disp.min = -2,
    disp.max = 2,
    slot = "scale.data",
    assay = "RNA",
    raster = TRUE,
    combine = TRUE
  ) + theme(legend.position='none')
print(plthm)
dev.off()



