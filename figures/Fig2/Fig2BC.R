library(Seurat)
library(dplyr)
library(ggplot2)
library(WriteXLS)
source("util.R")

genes <- c("Cre",
           "H2B-eGFP",
           "LacZ",
           "rtTA")
conditions <- c("Control",
                "Day12",
                "Day6",
                "Hour24",
                "Hour48_L",
                "Hour48_R",
                "Hour8"
)

object <- load_object("../save/scrna_phase_comparing.Rds")



mtx <- GetAssayData(object, slot="counts", assay="RNA")
dff <- data.frame(condition= character(0), Cre=integer(0), H2B_eGFP=integer(0), LacZ=integer(0), rtTA=integer(0))
for(cond in conditions){
  a_mtx <- mtx[, object$stage == cond]
  a_mtx[a_mtx > 0] <- 1
  adf <- data.frame(condition= cond, cells = ncol(a_mtx),
                                     Cre=sum(a_mtx["Cre", ]),
                                     H2B_eGFP=sum(a_mtx["H2B-eGFP", ]),
                                     LacZ=sum(a_mtx["LacZ", ]),
                                     rtTA=sum(a_mtx["rtTA", ]))
  dff <- rbind(dff, adf)
}


dff = dff %>% mutate(norm.Cre = dff$Cre * (dff$cells/sum(dff$cells)))
dff$norm.Cre <- round(prop.table(dff$norm.Cre) *100, 2)


dff = dff %>% mutate(norm.H2B_eGFP = dff$H2B_eGFP * (dff$cells/sum(dff$cells)))
dff$norm.H2B_eGFP <- round(prop.table(dff$norm.H2B_eGFP) *100, 2)

dff = dff %>% mutate(norm.LacZ = dff$LacZ * (dff$cells/sum(dff$cells)))
dff$norm.LacZ <- round(prop.table(dff$norm.LacZ) *100, 2)

dff = dff %>% mutate(norm.rtTA = dff$rtTA * (dff$cells/sum(dff$cells)))
dff$norm.rtTA <- round(prop.table(dff$norm.rtTA) *100, 2)

#dff <- dff %>% rename(cells=ncol)


WriteXLS(dff, "save/Markers_condition.xlsx", SheetNames="markers")

sdff <- dff[, c("condition", "norm.rtTA", "norm.LacZ", "norm.H2B_eGFP", "norm.rtTA")]
sdff <- reshape2::melt(sdff)

sdff$condition <- factor(sdff$condition, levels=c(
                                                  "Control",
                                                  "Hour8",
                                                  "Hour24",
                                                  "Hour48_L",
                                                  "Hour48_R",
                                                  "Day6",
                                                  "Day12")
)
pdf("save/Markers_condition.pdf")
ggplot(sdff, aes(x=condition, y=value, group=variable, fill=variable )) +
             geom_bar(position="dodge", stat="identity") +
             ylab("normalized proportion") +
             theme_minimal()
dev.off()

pdf("save/Markers_condition_proportion.pdf")

#sdff$condition <- factor(sdff$condition, levels=c("Control", "Hour24", "hour"))
ggplot(sdff, aes(x=variable, y=value, group=condition, fill=condition)) +
             geom_bar(stat="identity") +
             ylab("normalized proportion") +
             ggsci::scale_fill_simpsons() +
             ggsci::scale_color_simpsons() +
             theme_minimal()
dev.off()





