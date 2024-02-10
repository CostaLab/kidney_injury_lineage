library(Seurat)
library(ggplot2)

source("../util.R")

object <- load_object("../save/scrna_phase_comparing.Rds")


pdf("save/4weeks_umap.pdf")
DimPlot(object, group.by="removed_clusters", reduction="DEFAULT_UMAP") +
                          ggsci::scale_color_futurama() +
                          theme(legend.position='none') +
                          xlab("UMAP_1") +
                          ylab("UMAP_2")

DimPlot(object, group.by="removed_clusters", reduction="DEFAULT_UMAP") +
                          ggsci::scale_color_futurama() +
                          theme_void() +
                          theme(legend.position='none')+
                          ggtitle(NULL)

DimPlot(object, group.by="removed_clusters", label=T, reduction="DEFAULT_UMAP") + ggsci::scale_color_futurama()
dev.off()
