library(Seurat)

source("../../R/save_load_helper.R")


object <- load_object("./scrna_phase_comparing.Rds")

Idents(object) <- "name"


object <- subset(object, idents=c("week4_a", "week4_b"))
object$stage <- object$name


save_object(object, "scrna_rm_control.Rds")

