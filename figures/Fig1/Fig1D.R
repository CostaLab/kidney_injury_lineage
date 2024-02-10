library(Seurat)
library(ggplot2)
library(glue)
library(cowplot)
library(Rmagic)
library(WriteXLS)
library(patchwork)
library(Nebulosa)


genes <- c(
  "H2B-eGFP",
  "LacZ",
  "rtTA",
  "Cre"
)
source("../util.R")
scrna <- load_object("save/scrna_phase_comparing_subcluster_eGFP_STCs.Rds")

dir.create("figures")

pdf("figures/UMAPs.pdf")

p1 <- DimPlot(scrna, group.by="celltype",
               reduction="DEFAULT_UMAP",
               raster=T,
               label=F,
               cols=ggsci::pal_simpsons()(16))


p2 <- DimPlot(scrna, group.by="celltype2",
               reduction="DEFAULT_UMAP",
               raster=T,
               label=F,
               cols=ggsci::pal_simpsons()(16))


p3 <- DimPlot(scrna, group.by="annotation",
               reduction="DEFAULT_UMAP",
               raster=T,
               label=F,
               cols=ggsci::pal_simpsons()(16))

print(p1)
print(p1 + theme(legend.position='none'))
print(p2)
print(p2+ theme(legend.position='none'))
print(p3)
print(p3+ theme(legend.position='none'))

dev.off()



