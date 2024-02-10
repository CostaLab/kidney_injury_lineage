# The anaylsis of scRNA-seq data
We use single cell pipeline: https://github.com/CostaLab/scrna_seurat_pipeline

The analyses are according to the following configuration files:

1. config_1_B.R: first run of the pipeline
1. source_2_split.R: use resolutioin=0.5 but split sub-population into two clusters to 1&5 in resolution=0.8 to be x1 & x5.
1. config_3_B_split_0.5cluster2_Del57.R: remove 5 & 7 clusters
1. config_4_B_split_0.5cluster2_Del57_Del46x5.R: remove 4&6 and x5
1. source_5_rmControl.R: remove control cells
1. config_6_B_split_0.5cluster2_Del57_Del46x5_rmControl.R: re-run the analysis.
