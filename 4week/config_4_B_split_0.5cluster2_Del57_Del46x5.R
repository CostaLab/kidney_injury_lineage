### --------------Initail info----------------------------
PROJECT = "Acute Kidney Injury project 4weeks vs control" ## set project name
ORGAN = 'Kidney'           #For external annotation. Options: Blood, Heart, Intestine, Kidney
SPECIES = "Mouse"          #For external annotation. Options: Human, Mouse
MCA_NAME = "Kidney"        #For MCA annotation.      Options: check http://bis.zju.edu.cn/MCA/
HCL_NAME = "Kidney"        #For MCA annotation.      Options: check http://bis.zju.edu.cn/MCA/

# filtering params when create seurat object
MINCELLS  = 5
MINGENES  = 50


# filtering params when create seurat object
MINCELLS  = 5
MINGENES  = 50


INTEGRATION_OPTION = "seurat" ### or harmony

### -------------- Data SRC-----------------------------
ANNOTATION_EXTERNAL_FILE = "external/Human_and_mouse_cell_markers-Markers.tsv"

## If genesets you need are not included, please attach your geneset to the gmt.gz file.
MSigDB_GENESET_HUMAN_GMT_FILE  = "external/Human_msigdb.v7.2.symbols.gmt.gz"


data_src = c(
  "week4_a"      = "data/2022/cellranger/SID162686/outs/filtered_feature_bc_matrix",
  "week4_b"      = "data/2022/cellranger/SID162687/outs/filtered_feature_bc_matrix",
  "Control_a"    = "data/TX99_A006200089/cellranger/SID124771/outs/filtered_feature_bc_matrix",
  "Control_b"    = "data/TX139_A006850111/cellranger/SID142300/outs/filtered_feature_bc_matrix",
  "Control_c"    = "data/TX139_A006850111/cellranger/SID142301/outs/filtered_feature_bc_matrix"

)

##------------------ SET REPLICATE GROUP --------------
stage_lst = c(
        week4_a      =   "week4",
        week4_b      =   "week4",
        Control_a    =   "Control",
        Control_b    =   "Control",
        Control_c    =   "Control")



##------------------ Number of Cells for Doublet Detection --------------
# You have to set two variables.
# doublet_switch: This variable controls, if you want to detect doublets.
# Values:
# on      --> Detect doublets and remove them from the data set.
# off     --> Do not detect doublets.
# display --> Detect doublets, but do not remove them.
doublet_switch <- "on"

# This variable controls how you want to remove the doublets. If you set doublet_switch=off, you can ignore this variable.
# The percentage of doublets you expect depends on the number of cells you have. The more cells, the higher the percentage of doublets.
# 10X provides a table of the expected percentage of doublets and the number of loaded or recovered cells.
# The table is provided in static/DoubletEstimation10X.csv
# In this pipeline, you have the following options:
# 1) Use one set value, e.g. 0.05. This is the percentage of doublets, you want to remove from each sample.
# Example: doublet_lst <- 0.05

# 2) Use the number of cells in the Seurat object. The pipeline then uses the number of cells for each sample to estimate the proportion of doublets.
# Example: doublet_lst <- NULL

# 3) Set the percentage to be removed for each sample seperately.
# Example:
# doublet_lst = c(
#   A_MxCre = 0.05,
#   B_MxCre = 0.1,
#   C_Csnk  = 0.075,
#   D_Csnk  = 0.5
# )

# 4) Use the number of cells found by CellRanger for this sample.
# Then you have to provide the path to the CellRanger output file metrics_summary.csv.
# Example:
# doublet_lst = c(
#   A_MxCre = "/path/to/A_MxCre/metrics_summary.csv",
#   B_MxCre = "/path/to/B_MxCre/metrics_summary.csv",
#   C_Csnk  = "/path/to/C_Csnk/metrics_summary.csv",
#   D_Csnk  = "/path/to/D_Csnk/metrics_summary.csv"
# )

# 5) Set the number of cells, you want to use to calculate the proportion of doublets.
# Example:
# doublet_lst = c(
#   A_MxCre = 1024,
#   B_MxCre = 2048,
#   C_Csnk  = 4096,
#   D_Csnk  = 8192
# )

# 6) You can combine several options.
# Example:
# doublet_lst = c(
#   A_MxCre = "/path/to/A_MxCre/metrics_summary.csv", --> For this sample, we will use the number of cells after filtering by CellRanger.
#   B_MxCre = NULL, ------------------------------------> For this sample, we will use the number of cells in the current Seurat object.
#   C_Csnk  = 0.05, ------------------------------------> For this sample, we will simply remove 5% of the cells.
#   D_Csnk  = 1024 -------------------------------------> For this sample, we use this many cells to estimate the proportion of doublets.
# )



doublet_lst = c(
  "week4_a"      = "data/2022/cellranger/SID162686/outs/metrics_summary.csv",
  "week4_b"      = "data/2022/cellranger/SID162687/outs/metrics_summary.csv",
  "Control_a"    = "data/TX99_A006200089/cellranger/SID124771/outs/metrics_summary.csv",
  "Control_b"    = "data/TX139_A006850111/cellranger/SID142300/outs/metrics_summary.csv",
  "Control_c"    = "data/TX139_A006850111/cellranger/SID142301/outs/metrics_summary.csv"
)

## Phase_1, set 1 to regressout
preprocess_regressout = c("mito"       = 1,
                          "ribo"       = 0,
                          "cellcycle"  = 1
                         )


## Genesets candidate names please check external/MSigDB_names.txt
MSigDB_Geneset_names <- c(
    "NABA_COLLAGENS",
    "NABA_SECRETED_FACTORS",
    "NABA_ECM_GLYCOPROTEINS",
    "NABA_CORE_MATRISOME",
    "NABA_ECM_REGULATORS",
    "NABA_MATRISOME_ASSOCIATED",
    "NABA_ECM_AFFILIATED",
    "NABA_BASEMENT_MEMBRANES",
    "NABA_PROTEOGLYCANS",
    "NABA_MATRISOME"
)





#Analysis_phases
#1. scrna_phase_preprocess
#2. scrna_phase_clustering
#3. scrna_phase_comparing





### -------------- RUN PARAMETERS-----------------------------

## 0. omit,     1. calc & save,      2. load
conf = c(
       scrna_from_Del57           = 2,
       scrna_remove_clusters      = 1, ## remove clusters
       scrna_phase_comparing      = 1, ## DE GO pathway analysis etc. All rest calculating will be stored here
       scrna_cluster_annotation   = 0, ## Annotate clusters according to `cluster_annotation`
       scrna_clusterwise_xcell    = 0, ## remove cells of each cluster according distinct criterion
       scrna_del_mitogenes        = 0, ## !!!DANGEROUS, once deleted, never recovered!!!
       scrna_merge_clusters       = 0, ## merge clusters
       scrna_remove_recluster     = 0) ## remove clusters and recluster with default resolution

### ----------specific settings for some functions ----------------


## name[your operation name], value[dataframe which cluster, percentage to keep]
scrna_clusterwise_filtercell_settings <- list(
  "mito_cluster0,3,4,5"     =  data.frame(type="mito", max_pct=4, min_pct=0, cluster=c(0,3,4,5)),
  "ribo_cluster2,7_filter"  =  data.frame(type="ribo", max_pct=30, min_pct=0, cluster=c(2,7)),
  "mito_cluster6_filter"    =  data.frame(type="mito", max_pct=3, min_pct=0, cluster=6),
  "mito_cluster8_filter"    =  data.frame(type="mito", max_pct=5, min_pct=0, cluster=8)
)

## name[new cluster name], value[which clusters to merge together]
scrna_merge_clusters = list(
        "1+7" = c(1, 7),
        "2+6" = c(2, 6),
        "10+11+16" = c(10, 11, 16)
)

#clusters_split2

scrna_remove_clusters = c("4", "6", 'x5')
scrna_remove_recluster = c("4", "6", 'x5')


### cluster annotation
from_cluster_slot = "removed_clusters"
cluster_annotation <- c(
    "0" = "Vascular endothelial",
    "1" = "Fibroblasts 1",
    "2" = "Cardiomyocytes 1",
    "3" = "Endothelial Cells 1",
    "4" = "Macrophages",
    "6" = "Pericytes 1",
    "7" = "Cardiomyocytes 2",
    "8" = "Fibroblasts 2",
    "9" = "Cardiomyocytes 3",
    "10" = "Fibroblasts 3",
    "11" = "Lymphatic endothelial",
    "12" = "VSMCs",
    "13" = "Mesothelial cells",
    "14" = "Cardiomyocytes 4",
    "15" = "Lymphocytes",
    "16" = "T-cells 1",
    "18" = "Endothelial cells 2",
    "19" = "Pericytes 2",
    "20" = "T-cells 2",
    "21" = "Pericytes 2",
    "22" = "T-cells 3"
)

# Current options are the default option, you can change to your own
viz_conf = list(
  ## https://github.com/nanxstats/ggsci
  ##ggsci colorcode & availble colors:
                # "aaas":10 "d3":10 "futurama":12 "gsea":12 "igv":51
                # "jama":7 "jco":10 "lancet":9 "locuszoom":7 "material":10
                # "nejm":8 "npg":10 "rickandmorty":12 "simpsons":16 "startrek":7
                # "tron":7 "uchicago":9 "ucscgb":26

  cluster_color_option = "igv", ## ggsci, see above
  replicate_color_option = "simpsons", ## ggsci, see above
  neg_color = "#51C3CC",#colorBlindness::Blue2DarkOrange12Steps[2],
  pos_color = "#CC5800",#rev(colorBlindness::Blue2DarkOrange12Steps)[2],
  base_color = "lightgrey",#"lightgrey",
  neg_pos_divergent_palette = c('#1E8E99','#51C3CC','#99F9FF','#B2FCFF','#CCFEFF','#E5FFFF','#FFE5CC','#FFCA99','#FFAD65','#FF8E32','#CC5800','#993F00') #colorBlindness::Blue2DarkOrange12Steps
)


### MCA annotation Reference
#    "Arc-Me"                         "Bladder"
#    "Bone-Marrow"                    "Bone-Marrow_c-kit"
#    "Bone_Marrow_Mesenchyme"         "Brain"
#    "Bergman gliaBrain"              "E18-Brain"
#    "Embryonic-Mesenchyme"           "Embryonic-Stem-Cell"
#    "Female_Fetal_Gonad"             "Fetal_Brain"
#    "Fetal_Intestine"                "Fetal_Kidney"
#    "Fetal-liver"                    "Fetal_Lung"
#    "Fetal_Stomache"                 "Kidney"
#    "Liver"                          "Lung"
#    "Lung-Mesenchyme"                "Male_Fetal_Gonad"
#    "Mammary-Gland-Involution"       "Mammary-Gland-Lactation"
#    "Mammary-Gland-Pregrancy"        "Mammary-Gland-Virgin"
#    "Mesenchymal-Stem-Cell-Cultured" "Muscle"
#    "Neonatal_Brain"                 "Neonatal-Calvaria"
#    "Neonatal-Heart"                 "Neonatal-Muscle"
#    "Neonatal-Rib"                   "Neonatal-Skin"
#    "Ovary"                          "Pancreas"
#    "Peripheral_Blood"               "Placenta"
#    "Preimplantation-Embryo"         "Prostate"
#    "Retina"                         "Small-Intestine"
#    "Spleen"                         "Stomach"
#    "Testis"                         "Thymus"
#    "Trophoblast-Stem-Cell"          "Uterus"
