library(Seurat)
source("../../R/save_load_helper.R")


scrna <- load_object("./scrna_phase_comparing.Rds")

scrna$clusters_split2 <- as.character(scrna$seurat_clusters)

idx1 = which(scrna$integrated_snn_res.0.8 == "1")
idx5 = which(scrna$integrated_snn_res.0.8 == "5")


scrna$clusters_split2[idx1] = "x1"
scrna$clusters_split2[idx5] = "x5"
scrna$clusters_split2[scrna$clusters_split2=="2"] = 'x5' ## just simply set the rest 2 cluster2 cells to x5

save_object(scrna, "scrna_split.Rds", file_format="zstd")

