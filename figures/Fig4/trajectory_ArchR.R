mute <- suppressPackageStartupMessages
mute(library(Seurat))
mute(library(cowplot))
mute(library(dplyr))
mute(library(S4Vectors))
mute(library(SummarizedExperiment))
mute(library(gg3D))

addTrajectoryA <- function(
  scrna = NULL,
  name = "Trajectory",
  trajectory = c("9", "6", "5", "0", "2", "3", "4"),
  groupBy = "seurat_clusters",
  reducedDims = "diffusion",
  embedding = NULL,
  preFilterQuantile = 0.9,
  postFilterQuantile = 0.9,
  useAll = FALSE,
  dof = 250,
  spar = 1,
  force = FALSE,
  seed = 1
  ){

  if(!is.null(seed)) set.seed(seed)

  groupDF <- DataFrame(scrna@meta.data[, groupBy])
  rownames(groupDF) <- colnames(scrna)
  groupDF <- groupDF[groupDF[,1] %in% trajectory, , drop=F]


  if(sum(unique(groupDF[,1]) %in% trajectory)==0){
      stop
  }

  if(sum(unique(groupDF[,1]) %in% trajectory) < 3){
    stop
  }

  #mat <- scrna@reductions$diffusion@cell.embeddings
  mat <- scrna@reductions[[reducedDims]]@cell.embeddings
  mat <- mat[rownames(groupDF),,drop = FALSE]
  ######################################################
  #Filter Outliers
  ######################################################

  filterObj <- lapply(seq_along(trajectory), function(x){

      #Subset
      groupsx <- rownames(groupDF)[groupDF[,1]==trajectory[x]]
      matx <- mat[groupsx,,drop = FALSE]

      #Filter Distance
      matMeanx <- colMeans(matx)
      diffx <- sqrt(colSums((t(matx) - matMeanx)^2))
      idxKeep <- which(diffx <= quantile(diffx, preFilterQuantile))

      #Filter
      list(mat = matx[idxKeep,,drop=FALSE], groups = groupsx[idxKeep])

  })

  matFilter <- lapply(seq_along(filterObj), function(x) filterObj[[x]]$mat) %>% Reduce("rbind", .)
  groupsFilter <- groupDF[lapply(seq_along(filterObj), function(x) filterObj[[x]]$groups) %>% Reduce("c", .),,drop=FALSE]

  ######################################################
  #Now Initial Alignment
  ######################################################

  initialTime <- lapply(seq_along(trajectory), function(x){

      groupsx <- rownames(groupsFilter)[groupsFilter[,1] == trajectory[x]]
      matx <- matFilter[groupsx,,drop = FALSE]

      #Get Differences
      if(x != length(trajectory)){
          groupsxp1 <- rownames(groupsFilter)[groupsFilter[,1] == trajectory[x + 1]]
          meanx <- colMeans(matFilter[groupsxp1,,drop = FALSE])
          diffx <- sqrt(colSums((t(matx) - meanx)^2))
          timex <- (1 - ArchR:::.getQuantiles(diffx)) + x
      }else{
          groupsxm1 <- rownames(groupsFilter)[groupsFilter[,1] == trajectory[x - 1]]
          meanx <- colMeans(matFilter[groupsxm1,,drop = FALSE])
          diffx <- sqrt(colSums((t(matx) - meanx)^2))
          timex <- ArchR:::.getQuantiles(diffx) + x
      }

      timex

  }) %>% unlist

  ######################################################
  #Fit Cubic Splines
  ######################################################

  matSpline <- lapply(seq_len(ncol(matFilter)), function(x){
    tryCatch({
      stats::smooth.spline(
          x = initialTime,
          y = matFilter[names(initialTime), x],
          df = dof,
          spar = spar
      )[[2]]
    }, error = function(e){
      errorList <- list(
        it = x,
        x = initialTime,
        y = matFilter[names(initialTime), x],
        df = dof,
        spar = spar
      )

    })
  }) %>% Reduce("cbind",.) %>% data.frame()

  ######################################################
  # 1. KNN Fit vs Actual
  ######################################################

  knnObj <- nabor::knn(
      data =  matSpline,
      query = mat,
      k = 3
  )

  #Estimate place along trajectory
  knnIdx <- knnObj[[1]]
  knnDist <- knnObj[[2]]
  knnDiff <- ifelse(knnIdx[,2] > knnIdx[,3], 1, -1)
  knnDistQ <- ArchR:::.getQuantiles(knnDist[,1])

  #Filter Outlier Cells to Trajectory for High Resolution
  idxKeep <- which(knnDist[,1] <= quantile(knnDist[,1], postFilterQuantile))
  dfTrajectory <- DataFrame(
      row.names = rownames(mat),
      Distance = knnDist[, 1],
      DistanceIdx = knnIdx[, 1] + knnDiff * knnDistQ
  )[idxKeep, , drop = FALSE]


  dfTrajectory3 <- dfTrajectory

  dfTrajectory3$Trajectory <- 100 * ArchR:::.getQuantiles(dfTrajectory3[,2])
  nas <- rep(NA, dim(scrna)[2])
  names(nas) <- colnames(scrna)
  nas[rownames(dfTrajectory3)] <- dfTrajectory3$Trajectory
  scrna@meta.data[, name] <- nas
  return(scrna)
}

plotTrajectoryA <- function(
  scrna = NULL,
  trajectory = "Trajectory",
  embedding = "diffusion",
  name = "Trajectory",
  log2Norm = NULL,
  size = 0.2,
  rastr = TRUE,
  quantCut = c(0.01, 0.99),
  quantHex = 0.5,
  discreteSet = NULL,
  continuousSet = NULL,
  randomize = TRUE,
  keepAxis = FALSE,
  baseSize = 6,
  addArrow = TRUE,
  plotAs = NULL,
  smoothWindow = 5
  ){
  dfT <- DataFrame(scrna@meta.data[, trajectory])
  rownames(dfT) <- colnames(scrna)
  idxRemove <- which(is.na(dfT[,1]))
  df <- DataFrame(scrna@reductions[[embedding]]@cell.embeddings)
  dfT <- cbind(df, dfT[rownames(df),])
  colnames(dfT) <- c("DM1", "DM2", "PseudoTime")
  plotParams <- list()
  plotParams$x <- df[,1]
  plotParams$y <- df[,2]
  plotParams$title <- paste0(embedding, " of ", stringr::str_split(colnames(df)[1],pattern="#",simplify=TRUE)[,1])
  plotParams$baseSize <- baseSize
  plotParams$color <- as.vector(dfT$PseudoTime)
  plotParams$discrete <- ArchR:::.isDiscrete(plotParams$color)
  plotParams$continuousSet <- "horizonExtra"
  plotParams$discreteSet <- "stallion"
  plotParams$title <- paste(plotParams$title, " colored by\ncolData : ", name)

  .summarizeHex <- function(x = NULL){
    quantile(x, quantHex, na.rm = TRUE)
  }
  if(is.null(plotAs)){
      plotAs <- "hexplot"
  }
  plotParams$fun <- .summarizeHex
  plotParams$xlabel <- "DM1"
  plotParams$ylabel <- "DM2"

  if(!is.null(continuousSet)){
    plotParams$continuousSet <- continuousSet
  }
  if(!is.null(continuousSet)){
    plotParams$discreteSet <- discreteSet
  }
  plotParams$rastr <- rastr
  plotParams$size <- size
  plotParams$randomize <- randomize
  plotParams$color <- as.vector(plotParams$color)


  out <- do.call(ArchR::ggPoint, plotParams)
  out <- out + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
       dfT$value <- plotParams$color
  dfT <- dfT[order(dfT$PseudoTime), ]
  dfT <- dfT[!is.na(dfT$PseudoTime), ]

  colnames(dfT) <- c("x","y","PseudoTime","value")
  dfT <- as.data.frame(dfT)

  if(addArrow){
    dfArrow <-  split(dfT, floor(dfT$PseudoTime / 1.01)) %>%
                  lapply(colMeans) %>%
                  Reduce("rbind",.) %>%
                  data.frame

    dfArrow$x <- ArchR:::.centerRollMean(dfArrow$x, smoothWindow)
    dfArrow$y <- ArchR:::.centerRollMean(dfArrow$y, smoothWindow)
    dfArrow <- rbind(dfArrow, dfT[nrow(dfT), ,drop = FALSE])

    out <- out + geom_path(
            data = dfArrow, aes(x, y, color=NULL), size= 1,
            arrow = arrow(type = "open", length = unit(0.1, "inches"))
          )
  }

  return(out)
}

plotTrajectoryMulti <- function(
  scrna = NULL,
  trajectory = c("traj_1", "traj_3", "traj_4"),
  embedding = "diffusion",
  name = "Trajectory",
  log2Norm = NULL,
  size = 0.2,
  rastr = TRUE,
  quantCut = c(0.01, 0.99),
  quantHex = 0.5,
  discreteSet = NULL,
  continuousSet = NULL,
  randomize = TRUE,
  keepAxis = FALSE,
  baseSize = 6,
  addArrow = TRUE,
  plotAs = NULL,
  smoothWindow = 5
  ){

  dfT.list <- list()
  for(i in 1:length(trajectory)){
      dfT <- DataFrame(scrna@meta.data[, trajectory[i]])
      rownames(dfT) <- colnames(scrna)
      idxRemove <- which(is.na(dfT[,1]))
      df <- DataFrame(scrna@reductions[[embedding]]@cell.embeddings)
      dfT <- cbind(df, dfT[rownames(df),])
      colnames(dfT) <- c("DM1", "DM2", "PseudoTime")
      dfT.list[[i]] <- dfT

  }

  dfT <- DataFrame(scrna@meta.data[, trajectory[1]])
  rownames(dfT) <- colnames(scrna)
  idxRemove <- which(is.na(dfT[,1]))
  df <- DataFrame(scrna@reductions[[embedding]]@cell.embeddings)
  dfT <- cbind(df, dfT[rownames(df),])
  colnames(dfT) <- c("DM1", "DM2", "PseudoTime")
  #dfT.list[[i]] <- dfT
  #str(dfT)

  pseudoTime <- pseudo_vector_merge(lapply(dfT.list, function(x)x$PseudoTime))
  #str(pseudoTime)

  plotParams <- list()
  plotParams$x <- df[,1]
  plotParams$y <- df[,2]
  plotParams$title <- paste0(embedding, " of ", stringr::str_split(colnames(df)[1],pattern="#",simplify=TRUE)[,1])
  plotParams$baseSize <- baseSize
  #plotParams$color <- as.vector(dfT$PseudoTime)
  plotParams$color <- pseudoTime
  plotParams$discrete <- ArchR:::.isDiscrete(plotParams$color)
  plotParams$continuousSet <- "horizonExtra"
  plotParams$discreteSet <- "stallion"
  plotParams$title <- paste(plotParams$title, " colored by\ncolData : ", name)

  .summarizeHex <- function(x = NULL){
    quantile(x, quantHex, na.rm = TRUE)
  }
  if(is.null(plotAs)){
      plotAs <- "hexplot"
  }
  plotParams$fun <- .summarizeHex
  plotParams$xlabel <- "DM1"
  plotParams$ylabel <- "DM2"

  if(!is.null(continuousSet)){
    plotParams$continuousSet <- continuousSet
  }
  if(!is.null(continuousSet)){
    plotParams$discreteSet <- discreteSet
  }
  plotParams$rastr <- rastr
  plotParams$size <- size
  plotParams$randomize <- randomize
  plotParams$color <- as.vector(plotParams$color)


  out <- do.call(ArchR::ggPoint, plotParams)
  out <- out + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  for(i in 1:length(trajectory)){
    dfT <- dfT.list[[i]]
    dfT$value <- plotParams$color
    dfT <- dfT[order(dfT$PseudoTime), ]
    dfT <- dfT[!is.na(dfT$PseudoTime), ] ##

    colnames(dfT) <- c("x","y","PseudoTime","value")
    dfT <- as.data.frame(dfT)


    if(addArrow){
      dfArrow <-  split(dfT, floor(dfT$PseudoTime / 1.01)) %>%
                    lapply(colMeans) %>%
                    Reduce("rbind",.) %>%
                    data.frame

      dfArrow$x <- ArchR:::.centerRollMean(dfArrow$x, smoothWindow)
      dfArrow$y <- ArchR:::.centerRollMean(dfArrow$y, smoothWindow)
      dfArrow <- rbind(dfArrow, dfT[nrow(dfT), ,drop = FALSE])

      out <- out + geom_path(
              data = dfArrow, aes(x, y, color=NULL), size= 1,
              arrow = arrow(type = "open", length = unit(0.1, "inches"))
            )
    }

  }
  return(out)
}


plotTrajectory3D <- function(
  scrna = NULL,
  trajectory = "Trajectory",
  embedding = "diffusion",
  name = "Trajectory",
  log2Norm = NULL,
  size = 0.2,
  rastr = TRUE,
  quantCut = c(0.01, 0.99),
  quantHex = 0.5,
  discreteSet = NULL,
  continuousSet = NULL,
  randomize = TRUE,
  keepAxis = FALSE,
  baseSize = 6,
  addArrow = TRUE,
  plotAs = NULL,
  smoothWindow = 5,
  theta = 180,
  phi = 0
  ){
  dfT <- DataFrame(scrna@meta.data[, trajectory])
  rownames(dfT) <- colnames(scrna)
  idxRemove <- which(is.na(dfT[,1]))
  df <- DataFrame(scrna@reductions[[embedding]]@cell.embeddings)
  dfT <- cbind(df, dfT[rownames(df),])
  colnames(dfT) <- c("DM1", "DM2","DM3", "PseudoTime")
  plotParams <- list()
  plotParams$x <- df[,1]
  plotParams$y <- df[,2]
  plotParams$z <- df[,3]
  plotParams$title <- paste0(embedding, " of ", stringr::str_split(colnames(df)[1],pattern="#",simplify=TRUE)[,1])
  plotParams$baseSize <- baseSize
  plotParams$color <- as.vector(dfT$PseudoTime)
  plotParams$discrete <- ArchR:::.isDiscrete(plotParams$color)
  plotParams$continuousSet <- "horizonExtra"
  plotParams$discreteSet <- "stallion"
  plotParams$title <- paste(plotParams$title, " colored by\ncolData : ", name)

  .summarizeHex <- function(x = NULL){
    quantile(x, quantHex, na.rm = TRUE)
  }
  if(is.null(plotAs)){
      plotAs <- "hexplot"
  }
  plotParams$fun <- .summarizeHex
  plotParams$xlabel <- "DM1"
  plotParams$ylabel <- "DM2"
  plotParams$zlabel <- "DM3"

  if(!is.null(continuousSet)){
    plotParams$continuousSet <- continuousSet
  }
  if(!is.null(continuousSet)){
    plotParams$discreteSet <- discreteSet
  }
  plotParams$rastr <- rastr
  plotParams$size <- size
  plotParams$randomize <- randomize
  plotParams$color <- as.vector(plotParams$color)


  #out <- do.call(ArchR::ggPoint, plotParams)
  df3d <- data.frame(x=plotParams$x, y=plotParams$y, z=plotParams$z, color=plotParams$color)
  out <- ggplot(df3d, aes(x=x, y=y, z=z, color=color)) +    
          axes_3D(theta=theta, phi=phi) +
          stat_3D(theta=theta, phi=phi, geom="point",size=1, alpha = 3/5) +
          labs_3D(theta=theta, phi=phi, 
            labs=c("x", "y", "z"), 
            angle=c(0,0,0),
            hjust=c(0,2,2), 
            vjust=c(2,2,-2))


  out <- out + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  out <- out + theme_void() #+

  dfT$value <- plotParams$color
  dfT <- dfT[order(dfT$PseudoTime), ]
  dfT <- dfT[!is.na(dfT$PseudoTime), ]

  colnames(dfT) <- c("x","y","z","PseudoTime","value")
  dfT <- as.data.frame(dfT)

  if(addArrow){
    dfArrow <-  split(dfT, floor(dfT$PseudoTime / 1.01)) %>%
                  lapply(colMeans) %>%
                  Reduce("rbind",.) %>%
                  data.frame

    dfArrow$x <- ArchR:::.centerRollMean(dfArrow$x, smoothWindow)
    dfArrow$y <- ArchR:::.centerRollMean(dfArrow$y, smoothWindow)
    dfArrow$z <- ArchR:::.centerRollMean(dfArrow$z, smoothWindow)
    dfArrow <- rbind(dfArrow, dfT[nrow(dfT), ,drop = FALSE])

    out <- out + stat_3D(data=dfArrow,mapping=aes(x=x, y=y, z=z, color=NULL), 
                      theta=theta, phi=phi, geom="path",size=1,
                      arrow = arrow(type = "open", length = unit(0.1, "inches"))) 


  }

  return(out)
}




### haven't implemented because of the imputation
plotTrajectoryB <- function(
  scrna = NULL,
  trajectory = "Trajectory",
  embedding = "diffusion",
  name = "Cebpb", ## gene name
  log2Norm = NULL,
  size = 0.2,
  rastr = TRUE,
  quantCut = c(0.01, 0.99),
  quantHex = 0.5,
  discreteSet = NULL,
  continuousSet = NULL,
  randomize = TRUE,
  keepAxis = FALSE,
  baseSize = 6,
  addArrow = TRUE,
  plotAs = NULL,
  smoothWindow = 5
  ){
  dfT <- DataFrame(scrna@meta.data[, trajectory])
  rownames(dfT) <- colnames(scrna)
  idxRemove <- which(is.na(dfT[,1]))
  df <- DataFrame(scrna@reductions$diffusion@cell.embeddings)
  dfT <- cbind(df, dfT[rownames(df),])
  colnames(dfT) <- c("DM1", "DM2", "PseudoTime")
  plotParams <- list()
  plotParams$x <- df[,1]
  plotParams$y <- df[,2]
  plotParams$title <- paste0(embedding, " of ", stringr::str_split(colnames(df)[1],pattern="#",simplify=TRUE)[,1])
  plotParams$baseSize <- baseSize
  plotParams$color <- as.vector(dfT$PseudoTime)
  plotParams$discrete <- ArchR:::.isDiscrete(plotParams$color)
  plotParams$continuousSet <- "horizonExtra"
  plotParams$discreteSet <- "stallion"
  plotParams$title <- paste(plotParams$title, " colored by\ncolData : ", name)
  if(is.null(plotAs)){
      plotAs <- "hexplot"
  }
  plotParams$xlabel <- "DM1"
  plotParams$ylabel <- "DM2"

  if(!is.null(continuousSet)){
    plotParams$continuousSet <- continuousSet
  }
  if(!is.null(continuousSet)){
    plotParams$discreteSet <- discreteSet
  }
  plotParams$rastr <- rastr
  plotParams$size <- size
  plotParams$randomize <- randomize

  mat <- matrix(as.vector(plotParams$color), nrow = 1)
  colnames(mat) <- rownames(df)
  plotParams$color <- ArchR::imputeMatrix(mat = mat, imputeWeights = imputeWeights, logFile = "")

  plotParams$color <- ArchR:::.quantileCut(plotParams$color, min(quantCut), max(quantCut))


  dfT$value <- plotParams$color
  dfT <- dfT[order(dfT$PseudoTime), ]
  dfT <- dfT[!is.na(dfT$PseudoTime), ]

  colnames(dfT) <- c("x","y","PseudoTime","value")
  dfT <- as.data.frame(dfT)

  #Plot Pseudo-Time
  out2 <- ArchR:::ggPoint(
    x = dfT$PseudoTime,
    y = dfT$value,
    color = dfT$PseudoTime,
    discrete = FALSE,
    xlabel = "PseudoTime",
    ylabel = name,
    pal = plotParams$pal,
    ratioYX = 0.5,
    rastr = TRUE
  ) + geom_smooth(color = "black")

  attr(out2, "ratioYX") <- 0.5

  return(out2)
}



getTrajectoryA <- function(scrna = scrna,
  name = "Trajectory",
  useMatrix = "GeneScoreMatrix",
  groupEvery = 1,
  log2Norm = TRUE,
  scaleTo = 10000,
  smoothWindow = 11
){

    trajectory <- DataFrame(scrna@meta.data[, name])
    rownames(trajectory) <- colnames(scrna)
    trajectory <- trajectory[!is.na(trajectory[,1]),,drop=FALSE]

    breaks <- seq(0, 100, groupEvery)

    groupList <- lapply(seq_along(breaks), function(x){
          if(x == 1){
              NULL
          }else{
              rownames(trajectory)[which(trajectory[,1] > breaks[x - 1] & trajectory[,1] <= breaks[x])]
          }
    })[-1]
    names(groupList) <- paste0("T.", breaks[-length(breaks)], "_", breaks[-1])

    featureDF <-data.frame(seqnames = "z",
                             idx = seq_len(nrow(scrna)),
                             name = rownames(scrna),
                             stringsAsFactors = FALSE)


    seqnames <- unique(featureDF$seqnames)
    rownames(featureDF) <- paste0("f", seq_len(nrow(featureDF)))
    cellNames <- unlist(groupList, use.names = FALSE) ### UNIQUE here? doublet check QQQ
    matChr <- matrix(0, nrow = nrow(featureDF), ncol = length(groupList))
    colnames(matChr) <- names(groupList)
    rownames(matChr) <- rownames(featureDF)
    maty <- GetAssayData(scrna, slot="counts", assay="RNA")
    for(z in seq_along(groupList)){

      #Check Cells In Group
      cellsGroupz <- groupList[[z]]
      idx <- BiocGenerics::which(colnames(maty) %in% cellsGroupz)

      #If In Group RowSums
      if(length(idx) > 0){
        matChr[,z] <- matChr[,z] + Matrix::rowSums(maty[,idx,drop=FALSE])
      }
    }

    groupMat <- matChr

      if(!is.null(scaleTo)){
        if(any(groupMat < 0)){
          message("Some values are below 0, this could be a DeviationsMatrix in which scaleTo should be set = NULL.\nContinuing without depth normalization!")
        }else{
          groupMat <- t(t(groupMat) / colSums(groupMat)) * scaleTo
        }
      }

      if(log2Norm){
        if(any(groupMat < 0)){
          message("Some values are below 0, this could be a DeviationsMatrix in which log2Norm should be set = FALSE.\nContinuing without log2 normalization!")
        }else{
          groupMat <- log2(groupMat + 1)
        }
      }

      if(!is.null(smoothWindow)){

        message("Smoothing...")
        smoothGroupMat <- as.matrix(t(apply(groupMat, 1, function(x) ArchR:::.centerRollMean(x, k = smoothWindow))))
        colnames(smoothGroupMat) <- paste0(colnames(groupMat))
        colnames(groupMat) <- paste0(colnames(groupMat))

        #Create SE
        seTrajectory <- SummarizedExperiment(
            assays = SimpleList(
              smoothMat = as.matrix(smoothGroupMat),
              mat = as.matrix(groupMat)
            ),
            rowData = featureDF
        )
        if("name" %in% colnames(featureDF)){
          rownames(seTrajectory) <- paste0(featureDF$seqnames, ":", featureDF$name)
        }else{
          rownames(seTrajectory) <- paste0(featureDF$seqnames, ":", featureDF$start, "_", featureDF$end)
        }

      }else{

        colnames(groupMat) <- paste0(colnames(groupMat))

        #Create SE
        seTrajectory <- SummarizedExperiment(
            assays = SimpleList(
              mat = as.matrix(groupMat)
            ),
            rowData = featureDF
        )
        if("name" %in% colnames(featureDF)){
          rownames(seTrajectory) <- paste0(featureDF$seqnames, ":", featureDF$name)
        }else{
          rownames(seTrajectory) <- paste0(featureDF$seqnames, ":", featureDF$start, "_", featureDF$end)
        }
  }
      metadata(seTrajectory)$Params <- list(
        useMatrix = useMatrix,
        matrixClass = "Sparse.Assays.Matrix",
        scaleTo = scaleTo,
        log2Norm = log2Norm,
        smoothWindow = smoothWindow,
        date = Sys.Date()
      )
       seTrajectory
}
