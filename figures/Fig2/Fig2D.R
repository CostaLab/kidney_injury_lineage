library(Seurat)
library(ggplot2)
library(glue)
library(cowplot)
library(Rmagic)
library(WriteXLS)
library(patchwork)
library(Nebulosa)

#H2B-eGFP, LacZ, rt-TA and Cre

genes <- c(
  "H2B-eGFP",
  "LacZ",
  "rtTA",
  "Cre"
)
source("../util.R")
scrna <- load_object("save/scrna_phase_comparing.Rds")

dir.create("figures")


conds <- unique(scrna$stage)
pdf("figures/genes/featureplot.pdf")
DefaultAssay(scrna) <- "RNA"
p <- FeaturePlot(scrna, cols=c("lightgrey", "red"),
                        feature=genes,
                        reduction="DEFAULT_UMAP",
                        slot="data",
                        raster=F,
                        order=T) + plot_annotation(title="All Data")
print(p)

for (cond in conds){
  idx <- which(scrna$stage == cond)
  cells <- colnames(scrna)[idx]
  p <- FeaturePlot(scrna, cols=c("lightgrey", "red"),
                          cells = cells,
                          feature=genes,
                          reduction="DEFAULT_UMAP",
                          slot="data",
                          raster=F,
                          order=T) + plot_annotation(title=cond)
  print(p)

  #p <- plot_density(subset(scrna, cells=cells), features=genes) + plot_annotation(title=cond) + plot_layout(ncol=2)
  #print(p)
}

dev.off()



# conditions
cond_order<- c("Control", "Hour8", "Hour24", "Hour48_L", "Day6")
Idents(scrna) <- "stage"
cond_scrna <- subset(scrna, subset = stage %in% cond_order)
cond_scrna$stage <- factor(cond_scrna$stage, levels=cond_order)


pdf("figures/genes/Fig2DE_featureplot_genes.pdf", width=20, height=5)
DefaultAssay(cond_scrna) <- "RNA"

for(gene in genes){
  ps <- FeaturePlot(cond_scrna, cols=c("lightgrey", "red"),
                          feature=gene,
                          reduction="DEFAULT_UMAP",
                          split.by = "stage",
                          slot="data",
                          raster=F,
                          order=T,
                          combine = F)
  ps <- lapply(1:length(ps), function(i) ps[[i]] + theme_void() + ggtitle(cond_order[i]))
  p <- cowplot::plot_grid(plotlist=ps, nrow=1) + plot_annotation(title=gene)
  print(p)

  ps <- lapply(1:length(ps), function(i) ps[[i]] + theme(legend.position = "none"))
  p <- cowplot::plot_grid(plotlist=ps, nrow=1) + plot_annotation(title=gene)
  print(p)

}
dev.off()





scrna <- magic(scrna, genes=rownames(genes))
DefaultAssay(scrna) <- "MAGIC_RNA"
pdf("figures/Violinplot_magic_imputed.pdf")
for(g in genes){
  message(date(), " ", g)
  p1 <- VlnPlot(scrna, features=g, group.by="annotation", pt.size=0, cols=ggsci::pal_simpsons()(12))
  p2 <- VlnPlot(scrna, features=g, group.by="annotation", split.by="stage", pt.size=0, cols=ggsci::pal_simpsons()(12))
  p3 <- VlnPlot(scrna, features=g, group.by="stage", pt.size=0, cols=ggsci::pal_simpsons()(12))
  print(p1)
  print(p2)
  print(p3)
}
dev.off()


DefaultAssay(scrna) <- "RNA"
pdf("figures/Violinplot.pdf")
for(g in genes){
  message(date(), " ", g)
  p1 <- VlnPlot(scrna, features=g, group.by="annotation", pt.size=0, cols=ggsci::pal_simpsons()(12))
  p2 <- VlnPlot(scrna, features=g, group.by="annotation", split.by="stage", pt.size=0, cols=ggsci::pal_simpsons()(12))
  p3 <- VlnPlot(scrna, features=g, group.by="stage", pt.size=0, cols=ggsci::pal_simpsons()(12))
  print(p1)
  print(p2)
  print(p3)
}
dev.off()


df_list <- list()
for(g in genes){
 message(g)
 ev <- scrna@assays$RNA@counts[g, ] > 0
 dff <- data.frame(stage=scrna$stage, cellsExpressed=ev)
 df <- as.data.frame.matrix(table(dff$stage, dff$cellsExpressed))
 df_list[[g]] <- df
}

WriteXLS(df_list, "figures/counts.xlsx", row.names=T, SheetNames=names(df_list))


df_list <- list()
scrna$stage_clusters <- paste0(scrna$stage, "-", scrna$annotation)
for(g in genes){
 message(g)
 ev <- scrna@assays$RNA@counts[g, ] > 0
 dff <- data.frame(stage_clusters=scrna$stage_clusters, cellsExpressed=ev)
 df <- as.data.frame.matrix(table(dff$stage_clusters, dff$cellsExpressed))
 df_list[[g]] <- df
}
WriteXLS(df_list, "figures/counts_clusters.xlsx", row.names=T, SheetNames=names(df_list))


