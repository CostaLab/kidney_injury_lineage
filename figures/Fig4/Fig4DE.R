library(Seurat)
library(glue)
library(Nebulosa)
library(tidyverse)




object <- load_object("../save/scrna_phase_comparing.Rds")

DefaultAssay(object) <- "MAGIC_RNA"

feas <- c("Cd74", "Mif", "Mcat")


pdf("save/fig4DE.pdf")
for(fea in feas){

  p <- plot_density(object, features=fea, reduction="DEFAULT_UMAP")
  print(p)

}
dev.off()

