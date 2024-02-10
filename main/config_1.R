### --------------Initail info----------------------------
PROJECT = "Acute Kidney Injury project 3 replicate each condition" ## set project name
ORGAN = 'Kidney'           #For external annotation. Options: Blood, Heart, Intestine, Kidney
SPECIES = "Mouse"          #For external annotation. Options: Human, Mouse
MCA_NAME = "Kidney"        #For MCA annotation.      Options: check http://bis.zju.edu.cn/MCA/
HCL_NAME = "Kidney"        #For MCA annotation.      Options: check http://bis.zju.edu.cn/MCA/

# filtering params when create seurat object
MINCELLS  = 5 
MINGENES  = 50

### -------------- Data SRC-----------------------------
ANNOTATION_EXTERNAL_FILE = "external/Human_and_mouse_cell_markers-Markers.tsv" 

INTEGRATION_OPTION = "seurat" ### or harmon

#SID142302: 12 days after injury
#SID142303: 6 days after injury
#SID142307: 8 hours after injury
#SID142304: 24 hours after injury
#SID142305: right kidney, 48 h after injury
#SID142306: left kidney, 48h after injury
#SID142300, SID142301: controls
#SID126683: Hour48_L
#SID126684: Hour48_R


data_src = c( 
  "Day12_a"      = "data/TX92_scRNAv3_A006200068/cellranger/SID118110/outs/filtered_feature_bc_matrix",
  "Day12_b"      = "data/TX92_scRNAv3_A006200068/cellranger/SID118111/outs/filtered_feature_bc_matrix",
  "Day12_c"      = "data/TX139_A006850111/cellranger/SID142302/outs/filtered_feature_bc_matrix",
  "Day6_a"       = "data/TX92_scRNAv3_A006200068/cellranger/SID118112/outs/filtered_feature_bc_matrix",
  "Day6_b"       = "data/TX92_scRNAv3_A006200068/cellranger/SID118113/outs/filtered_feature_bc_matrix",
  "Day6_c"       = "data/TX139_A006850111/cellranger/SID142303/outs/filtered_feature_bc_matrix",
  "Hour48_R_a"   = "data/TX99_A006200089/cellranger/SID124769/outs/filtered_feature_bc_matrix",
  "Hour48_R_b"   = "data/TX139_A006850111/cellranger/SID142305/outs/filtered_feature_bc_matrix",
  "Hour48_R_c"   = "data/TX106_A006850073_ES/SID126684/outs/filtered_feature_bc_matrix",
  "Hour48_L_a"   = "data/TX99_A006200089/cellranger/SID124770/outs/filtered_feature_bc_matrix",
  "Hour48_L_b"   = "data/TX139_A006850111/cellranger/SID142306/outs/filtered_feature_bc_matrix",
  "Hour48_L_c"   = "data/TX106_A006850073_ES/SID126683/outs/filtered_feature_bc_matrix",
  "Hour24_a"     = "data/TX92_scRNAv3_A006200068/cellranger/SID118115/outs/filtered_feature_bc_matrix",
  "Hour24_b"     = "data/TX99_A006200089/cellranger/SID124773/outs/filtered_feature_bc_matrix", 
  "Hour24_c"     = "data/TX139_A006850111/cellranger/SID142304/outs/filtered_feature_bc_matrix",
  "Hour8_a"      = "data/TX92_scRNAv3_A006200068/cellranger/SID118114/outs/filtered_feature_bc_matrix",
  "Hour8_b"      = "data/TX92_scRNAv3_A006200068/cellranger/SID118116/outs/filtered_feature_bc_matrix",
  "Hour8_c"      = "data/TX139_A006850111/cellranger/SID142307/outs/filtered_feature_bc_matrix",
  "Control_a"    = "data/TX99_A006200089/cellranger/SID124771/outs/filtered_feature_bc_matrix",
  "Control_b"    = "data/TX139_A006850111/cellranger/SID142300/outs/filtered_feature_bc_matrix",
  "Control_c"    = "data/TX139_A006850111/cellranger/SID142301/outs/filtered_feature_bc_matrix"

)

#A_MxCre B_MxCre  C_Csnk  D_Csnk 
##------------------ SET REPLICATE GROUP --------------
stage_lst = c(
        Day12_a      =   "Day12",
        Day12_b      =   "Day12",
        Day12_c      =   "Day12",
        Day6_a       =   "Day6",
        Day6_b       =   "Day6",
        Day6_c       =   "Day6",
        Hour48_R_a   =   "Hour48_R",
        Hour48_R_b   =   "Hour48_R",
        Hour48_R_c   =   "Hour48_R",
        Hour48_L_a   =   "Hour48_L",
        Hour48_L_b   =   "Hour48_L",
        Hour48_L_c   =   "Hour48_L",
        Hour24_a     =   "Hour24",
        Hour24_b     =   "Hour24",
        Hour24_c     =   "Hour24",
        Hour8_a      =   "Hour8",
        Hour8_b      =   "Hour8",
        Hour8_c      =   "Hour8",
        Control_a    =   "Control",
        Control_b    =   "Control",
        Control_c    =   "Control")




preprocess_regressout = c("mito"      = 1,
                          "ribo"      = 0,
                          "cellcycle" = 1)



#control, 8h, 24h,48hL, 6d, 12d, 48hR
stage_order <- c(
        "Control",
        "Hour8",
        "Hour24",
        "Hour48_L",
        "Day6",
        "Day12",
        "Hour48_R")

sample_order <- c(
        "Control_a",
        "Control_b",
        "Control_c",
        "Hour8_a",
        "Hour8_b",
        "Hour8_c",
        "Hour24_a",
        "Hour24_b",
        "Hour24_c",
        "Hour48_L_a",
        "Hour48_L_b",
        "Hour48_L_c",
        "Day6_a",
        "Day6_b",
        "Day6_c",
        "Day12_a",
        "Day12_b",
        "Day12_c",
        "Hour48_R_a",
        "Hour48_R_b",
        "Hour48_R_c")
### -------------- RUN PARAMETERS-----------------------------


### -------------- RUN PARAMETERS-----------------------------

## 0. omit,     1. calc & save,      2. load
conf = c(
       scrna_phase_preprocess     = 1, ## quality check and preprocessing before integration
       scrna_phase_clustering     = 1, ## integration & clustering
       scrna_phase_comparing      = 1, ## DE GO pathway analysis etc. All rest calculating will be stored here
       scrna_progeny              = 1,
       scrna_progeny_stage        = 1,
       scrna_cluster_annotation   = 0, ## Annotate clusters according to `cluster_annotation`
       scrna_clusterwise_xcell    = 0, ## remove cells of each cluster according distinct criterion
       scrna_del_mitogenes        = 0, ## !!!DANGEROUS, once deleted, never recovered!!!
       scrna_merge_clusters       = 0, ## merge clusters
       scrna_remove_clusters      = 0, ## remove clusters
       scrna_remove_recluster     = 0) ## remove clusters and recluster with default resolution



### ----------specific settings for some functions ----------------

scrna_merge_clusters = list(
        "1+7" = c(1, 7),
        "2+6" = c(2, 6),
        "10+11+16" = c(10, 11, 16)
)


scrna_remove_clusters = c(1, 3, 6)
scrna_remove_recluster = c(1, 3, 6)


# Current options are the default option, you can change to your own
viz_conf = list(
  ## see: https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html#the-color-scales
  # "magma" (or "A"), "inferno" (or "B"), "plasma" (or "C"), "viridis" (or "D") and "cividis" (or "E")
  cluster_color_option = "D", ## R viridis, see above
  replicate_color_option = "C", ## R viridis, see above
  neg_color = "#51C3CC",#colorBlindness::Blue2DarkOrange12Steps[2],
  pos_color = "#CC5800",#rev(colorBlindness::Blue2DarkOrange12Steps)[2],
  base_color = "lightgrey",#"lightgrey",
  neg_pos_divergent_palette = c('#1E8E99','#51C3CC','#99F9FF','#B2FCFF','#CCFEFF','#E5FFFF','#FFE5CC','#FFCA99','#FFAD65','#FF8E32','#CC5800','#993F00') #colorBlindness::Blue2DarkOrange12Steps
)

