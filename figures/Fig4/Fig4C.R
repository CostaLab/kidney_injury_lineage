mute <- suppressPackageStartupMessages
mute(library(destiny)) # diffusion maplibrary(Seurat) # Seurat3
mute(library(Seurat)) # Seurat3
mute(library(ArchR)) # For Heatmap
mute(library(parallelDist))
mute(library(schex))
mute(library(dplyr))
mute(library(reshape2))
mute(library(gg3D))
mute(library(plotly))


source("trajectory_ArchR.R")
source("../util.R")

## scrna is a Seurat 3 object
scrna <- load_object(file="../save/scrna_phase_comparing.Rds")

Idents(scrna) <- "annotation"

all_clusters <- droplevels(unique(scrna$annotation))
keep <- setdiff(all_clusters, c("PT4", "PT6", "PT7", "PT8"))

a_sub <- subset(scrna, idents=keep)


dm <- destiny::DiffusionMap(a_sub@reductions$INTE_PCA@cell.embeddings, verbose = TRUE)
embedding <- as.data.frame(dm@eigenvectors)[, c("DC1", "DC2")]
colnames(embedding) <- paste0("DM_", 1:2)
rownames(embedding) <- colnames(a_sub)
a_sub[["diffusion"]] <- CreateDimReducObject(embeddings = as.matrix(embedding), key = "DM_", assay = DefaultAssay(a_sub))

a_sub <- FindNeighbors(a_sub, reduction = "diffusion", dims = 1:2, verbose = T)
a_sub <- FindClusters(a_sub, resolution = 0.8, verbose = T) ##
a_sub$diffusion_clusters <- a_sub$seurat_clusters



a_sub <- make_hexbin(a_sub, nbins = 100,
    dimension_reduction = "diffusion", use_dims=c(1,2))

dir.create("viz")
a_sub <- FindNeighbors(a_sub, reduction = "diffusion", dims = 1:2, verbose = T)
a_sub <- FindClusters(a_sub, resolution = 0.8, verbose = T) ##
a_sub$diffusion_clusters <- a_sub$seurat_clusters
pdf("viz/diffusion2d_embeddings.pdf", width=14, height=10)

a_sub <- make_hexbin(a_sub, nbins = 100,
    dimension_reduction = "diffusion", use_dims=c(1,2))

DimPlot(a_sub, reduction="diffusion", pt.size=0.1, group.by="stage", label=T,cols=ggsci::pal_igv()(15))
plot_hexbin_meta(a_sub, col="stage", action="majority") + ggsci::scale_fill_igv()
DimPlot(a_sub, reduction="diffusion", pt.size=0.1, group.by="celltype", label=T,cols=ggsci::pal_igv()(15))
plot_hexbin_meta(a_sub, col="celltype", action="majority") + ggsci::scale_fill_igv()
DimPlot(a_sub, reduction="diffusion", pt.size=0.1, group.by="celltype2", label=T,cols=ggsci::pal_igv()(15))
plot_hexbin_meta(a_sub, col="celltype2", action="majority") + ggsci::scale_fill_igv()
DimPlot(a_sub, reduction="diffusion", pt.size=0.1, group.by="diffusion_clusters",label=T, cols=ggsci::pal_igv()(50))
plot_hexbin_meta(a_sub, col="diffusion_clusters", action="majority")+ ggsci::scale_fill_igv()

dev.off()


traj_1 <- c("42",
            "43",
            "41",
            "38",
            "33",
            "22",
            "28",
            "35",
            "31"
)


traj_2 <- c("17",
            "26",
            "15",
            "32",
            "7",
            "16",
            "4",
            "20",
            "34",
            "10",
            "0",
            "9",
            "31"
)



traj_3 <- c("42",
            "43",
            "41",
            "38",
            "37",
            "2",
            "8",
            "14",
            "18",
            "11",
            "17"
)


traj_4 <- c("42",
            "43",
            "41",
            "38",
            "37",
            "24",
            "25",
            "32",
            "15"
)


pdf("viz/trajectorys2d.pdf")

## clustering
# add trajectory based along a group of clusters

plotTrajectoryMulti(a_sub, name="traj_2", embedding="diffusion", trajectory=c("traj_1", "traj_3", "traj_4"))

a_sub <-  addTrajectoryA(a_sub, name="traj_1", trajectory = traj_1, groupBy="diffusion_clusters" )
p <- plotTrajectoryA(a_sub, name="traj_1", embedding="diffusion", trajectory="traj_1")
print(p)

a_sub <-  addTrajectoryA(a_sub, name="traj_2", trajectory = traj_2, groupBy="diffusion_clusters" )
p <- plotTrajectoryA(a_sub, name="traj_2", embedding="diffusion", trajectory="traj_2")
print(p)

a_sub <-  addTrajectoryA(a_sub, name="traj_3", trajectory = traj_3, groupBy="diffusion_clusters" )
p <- plotTrajectoryA(a_sub, name="traj_3", embedding="diffusion", trajectory="traj_3")
print(p)

a_sub <-  addTrajectoryA(a_sub, name="traj_4", trajectory = traj_4, groupBy="diffusion_clusters" )
p <- plotTrajectoryA(a_sub, name="traj_4", embedding="diffusion", trajectory="traj_4")
print(p)
dev.off()






pdf("viz/trajectory2d_heatmaps.pdf", height=14, width=10)
## plot trajectory heatmap
se <- getTrajectoryA(a_sub, name="traj_1")
ht <- ArchR::plotTrajectoryHeatmap(se, labelTop=100)
mtx <-  ArchR::plotTrajectoryHeatmap(se, labelTop=100, returnMatrix=T, maxFeatures=100)
genes <- rownames(mtx)
draw(ht, column_title = "traj_1")

## plot trajectory heatmap
se <- getTrajectoryA(a_sub, name="traj_2")
ht <- ArchR::plotTrajectoryHeatmap(se, labelTop=100)
mtx <-  ArchR::plotTrajectoryHeatmap(se, labelTop=100, returnMatrix=T, maxFeatures=100)
genes <- union(genes, rownames(mtx))
draw(ht, column_title = "traj_2")

## plot trajectory heatmap
se <- getTrajectoryA(a_sub, name="traj_3")
ht <- ArchR::plotTrajectoryHeatmap(se, labelTop=100)
mtx <-  ArchR::plotTrajectoryHeatmap(se, labelTop=100, returnMatrix=T, maxFeatures=100)
genes <- union(genes, rownames(mtx))
draw(ht, column_title = "traj_3")

## plot trajectory heatmap
se <- getTrajectoryA(a_sub, name="traj_4")
ht <- ArchR::plotTrajectoryHeatmap(se, labelTop=100)
mtx <-  ArchR::plotTrajectoryHeatmap(se, labelTop=100, returnMatrix=T, maxFeatures=100)
genes <- union(genes, rownames(mtx))
draw(ht, column_title = "traj_4")
dev.off()


DefaultAssay(a_sub) <- "MAGIC_RNA"

gene_mtx <- GetAssayData(a_sub, slot="data")
fgenes <- genes %>% stringr::str_sub(start=3)


pdf("viz/linePlots2d.pdf", width=10, height=5)
for(gene in fgenes){
  gene_vec <- gene_mtx[gene, ]
  dfT1 <- as.data.frame(a_sub$traj_1)
  colnames(dfT1) <- "traj_1"

  dfT1 <- dfT1 %>% dplyr::filter(!is.na(traj_1))
  dfT1 <- cbind(dfT1, data.frame(gene_vec[rownames(dfT1)]))
  colnames(dfT1) <- c("PseudoTime", "genescore")
  dfT1$name <- "traj_1"

  dfT3 <- as.data.frame(a_sub$traj_3)
  colnames(dfT3) <- "traj_3"
  dfT3 <- dfT3 %>% dplyr::filter(!is.na(traj_3))
  dfT3 <- cbind(dfT3, data.frame(gene_vec[rownames(dfT3)]))
  colnames(dfT3) <- c("PseudoTime", "genescore")
  dfT3$name <- "traj_3"

  dfT4 <- as.data.frame(a_sub$traj_4)
  colnames(dfT4) <- "traj_4"
  dfT4 <- dfT4 %>% dplyr::filter(!is.na(traj_4))
  dfT4 <- cbind(dfT4, data.frame(gene_vec[rownames(dfT4)]))
  colnames(dfT4) <- c("PseudoTime", "genescore")
  dfT4$name <- "traj_4"

 
  p <- ggplot() +
              geom_smooth(data=dfT1, aes(x=PseudoTime, y=genescore, color=name), method = "loess") +
              #geom_smooth(data=dfT2, aes(x=PseudoTime, y=genescore, color=name), method = "loess") +
              geom_smooth(data=dfT3, aes(x=PseudoTime, y=genescore, color=name), method = "loess") +
              geom_smooth(data=dfT4, aes(x=PseudoTime, y=genescore, color=name), method = "loess") +
              theme_cowplot()+
              labs(color="Legend") +
              ggtitle(gene)

  print(p)
}

dev.off()


saveRDS(a_sub, "save/trajectory_added.Rds")




