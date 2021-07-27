# Script to determine stimulation NN time

### Author: Vinay Kartha
### Contact: <vinay_kartha@g.harvard.edu>
### Affiliation: Buenrostro Lab, Department of Stem Cell and Regenerative Biology, Harvard University


library(dplyr)
library(chromVAR)
library(BuenColors)
library(ggplot2)


# Function to make dotPlot of nn stim time estimates, using Control 1h, and the stim 1h 6h and GI conditions
nnDotPlot <- function(df, # Pseudotime data frame (melted). Must have columns "Condition", "cellType", "variable", and "value"
                      stimType="IFN",
                      cellType=NULL){
  
  stopifnot(stimType %in% df$variable)
  #For this stim
  myConditions <- c("Control_1h",paste0(stimType,c("_1h","_6h","GolgiPlug_6h")))
  cat("Keeping conditions: ",myConditions,"\n")
  stopifnot(all(myConditions %in% df$Condition))
  
  if(is.null(cellType)){
    myCellType <- unique(df$cellType)
  } else {
    # Filter for cellType
    stopifnot(all(cellType %in% unique(df$cellType)))
    myCellType <- cellType
  }
  # Customize x-axis labels (grouping)
  myLabels <- c("1h","1h","6h","GI")
  names(myLabels) <- myConditions
  
  # Edited below line of code to plot scaled values instead of raw nn stim time estimates
  df %>% filter(Condition %in% myConditions & variable==paste0(stimType,".scaled") & cellType %in% myCellType) %>%
    ggplot( aes(x=factor(Condition),y=value,color=Condition)) + 
    geom_quasirandom(groupOnX = TRUE,size=0.8,shape=16) +
    # Add median value
    stat_summary(fun = "median", size= 0.3, geom = "crossbar",color="black",width=0.5)+
    #geom_boxplot(outlier.size = 0.1,outlier.color = NA,width=0.5,color="black",fill=NA,alpha=0.5)+
    theme_classic() + scale_y_continuous(expand=c(0,0)) +
    labs(y="Stim response time",x="Condition") + 
    scale_color_manual(values = myCols) + 
    facet_grid(~cellType) + 
    scale_x_discrete(labels=myLabels)+
    theme(legend.position = "none",
          #axis.text.x = element_text(angle=60,vjust=1,hjust = 1),
          axis.text = element_text(color="black"),
          strip.background = element_blank(),strip.text = element_text(size=11))
  
}


# Function to compute NN stim time using either ATAC/RNA dorcs
computeNNtime <- function(meta,# ATAC or RNA meta data (must have valid Condition column)
                          mat, # LSI matrix or smoothed normed DORC score matrix
                          isDORC=FALSE, # Is the above mat a DORC matrix? Assumes No (i.e. it is LSI/PCA).
                          K=50 # K to use to find NN stim time based on nearest neighbor labels
){
  
  # Do knn of unntreated cells to control+stim cells (leaving out Golgi here)
  untreatedCells <- meta$Condition %in% "Untreated"
  cellsToKeep <- !untreatedCells & !grepl("Golgi",meta$Condition)  # Includes controls, excludes golgis, untreated
  cat("Cells being considered for NN mapping\n")
  print(table(meta$Condition,cellsToKeep))
  
  # Conditions of all stim cells (minus Golgi and Untreated)
  stimConditions <- factor(meta$Condition[cellsToKeep])
  levels(stimConditions)
  
  # Give 0 weight to control cells, and 1 and 6 units to 1hr and 6hr cells in each stim, respectively
  stimWeights <- c("Control_1h"=0,"Control_6h"=0,
                   "IFN_1h"=1,"IFN_6h"=2, # We tried both 2/6 hr weight here
                   "LPS_1h"=1,"LPS_6h"=2, # We tried both 2/6 hr weight here
                   "PMA_1h"=1,"PMA_6h"=2) # We tried both 2/6 hr weight here
  
  stopifnot(all(names(stimWeights) %in% levels(stimConditions)))
  
  # Map knns for each query cell to each control/stim based on chromatin
  
  # If using sig smoothed DORC PC scores
  # PCA on DORC scores Takes a long time to run if not cached already
  # Should be cached already
  if(isDORC){
    cat("Determining PCA on DORC scores ..\n")
    pc.DORC <- cachePCA(dataSet=t(scale(mat)),
                        scale=FALSE,
                        center=FALSE,
                        cachePath = "./PCAcache/")
    
    
    DORCpcScores <- pc.DORC$x[,1:20]
    
    cat("PC score dimensions\n")
    print(dim(DORCpcScores))
    
    cat("Getting ",K," nearest neighbors for NN stim time ..\n")
    # All stim cells to all stim cells (minus golgi) map based on DORC PCs
    knn <- FNN::get.knnx(query = DORCpcScores,
                         data = DORCpcScores[cellsToKeep,],k = K)$nn.index  # Excludes  golgis and untreated from neighbor mapping (in ref)
  } else {
    
    cat("Getting ",K," nearest neighbors for NN stim time ..\n")
    # All stim cells to all stim cells (minus golgi) map based on LSI PCs
    knn <- FNN::get.knnx(query = mat,
                         data = mat[cellsToKeep,],k = K)$nn.index  # Excludes  golgis and untreated from neighbor mapping (in ref)
  }
  
  dim(knn)
  
  # For each cell, compute the mean stim response pseudotime based on the matched control/stim cells
  # Do this for each stim treatment condition
  
  pseudoTimelist <- list()
  for(i in 1:nrow(knn)){
    pseudoTimelist[[i]] <- table(stimConditions[knn[i,]]) * stimWeights/ncol(knn)
    
  }
  
  pseudoTime <- do.call('rbind',pseudoTimelist)
  colnames(pseudoTime)
  
  # Weighted sum of 1hr and 6hr for each stim condition
  # Each will be vector of size ncol(ATAC)
  IFN.pseudo <- rowSums(pseudoTime[,c("IFN_1h","IFN_6h")])
  LPS.pseudo <- rowSums(pseudoTime[,c("LPS_1h","LPS_6h")])
  PMA.pseudo <- rowSums(pseudoTime[,c("PMA_1h","PMA_6h")])
  
  pseudo.d <- data.frame("IFN"=IFN.pseudo,
                         "LPS"=LPS.pseudo,
                         "PMA"=PMA.pseudo,
                         Condition = as.character(meta$Condition), # This lets you filter cells by Condition
                         stringsAsFactors = FALSE)
  
  rownames(pseudo.d) <- rownames(meta)
  
  return(pseudo.d)
}


setwd("<data_analysis_folder>")

source("./code/misc_helper_stim.R")

# ATAC based NN time

# Load ATAC SE (paired)
# NOTE that we load the paired dataset here
ATAC.SE <- readRDS("./data/SE/ATAC_stim_paired.rds")
ATAC.meta <- as.data.frame(colData(ATAC.SE),stringsAsFactors=FALSE)


# V2 below is for weights 0-2; 
pseudoFileATAC <- "./data/nnTime/LSI_NN_stim_pseudotimev2.tsv"

if(file.exists(pseudoFileATAC)){
  pseudo.d <- read.table(pseudoFileATAC,sep="\t",header=TRUE,stringsAsFactors = FALSE) 
} else {

  # Load LSI on full ATAC data (not just paired)
  lsi <- read.table("./data/ArchR/scATAC_stim_LSI.tsv",sep="\t",header=TRUE,stringsAsFactors=FALSE,row.names = 1)
  head(lsi)
  
  stopifnot(all(rownames(ATAC.meta) %in% rownames(lsi)))
  
  # Subset to paired ATAC
  lsi <- lsi[colnames(ATAC.SE),]
  
  pseudo.d <- computeNNtime(meta = ATAC.meta,mat = lsi,K = 50) # Computed across all cells using condition labels
  
  # Re-scaled to 0-1
  pseudo.d$IFN.scaled <- scales::rescale(pseudo.d$IFN,to=c(0,1))
  pseudo.d$LPS.scaled <- scales::rescale(pseudo.d$LPS,to=c(0,1))
  pseudo.d$PMA.scaled <- scales::rescale(pseudo.d$PMA,to=c(0,1))

  stopifnot(identical(rownames(ATAC.meta),rownames(pseudo.d)))
  
  pseudo.d$cellType <- ATAC.meta$pairedLabel2
  write.table(pseudo.d,pseudoFileATAC,sep="\t",quote=FALSE)
}

stopifnot(identical(rownames(ATAC.meta),rownames(pseudo.d)))


# UMAP of NN stim time for each stim (Fig 4 in paper)
library(ggplot2)
library(BuenColors)

IFN_UMAP_time <- plotMarker2D(ATAC.meta[,c("UMAP1","UMAP2")],markerMat = t(pseudo.d[,paste0(c("IFN","LPS","PMA"),".scaled")]),markers = "IFN.scaled",labels = "IFN_NN_time",pointSize = 0.1,colorPalette = "brewer_heat",rasteRize = TRUE) + scale_color_gradientn(breaks=scales::pretty_breaks(n=1),colours = jdb_palette("brewer_heat"))
LPS_UMAP_time <- plotMarker2D(ATAC.meta[,c("UMAP1","UMAP2")],markerMat = t(pseudo.d[,paste0(c("IFN","LPS","PMA"),".scaled")]),markers = "LPS.scaled",labels = "LPS_NN_time",plotTitle = "LPS NN time",pointSize = 0.1,colorPalette = "brewer_blue",rasteRize = TRUE) + scale_color_gradientn(breaks=scales::pretty_breaks(n=1),colours = jdb_palette("brewer_blue"))
PMA_UMAP_time <- plotMarker2D(ATAC.meta[,c("UMAP1","UMAP2")],markerMat = t(pseudo.d[,paste0(c("IFN","LPS","PMA"),".scaled")]),markers = "PMA.scaled",labels = "PMA_NN_time",plotTitle = "PMA NN time",pointSize = 0.1,colorPalette = "brewer_green",rasteRize = TRUE) + scale_color_gradientn(breaks=scales::pretty_breaks(n=1),colours = jdb_palette("brewer_green"))

cowplot::plot_grid(IFN_UMAP_time,LPS_UMAP_time,PMA_UMAP_time,nrow=1,align="hv")

d.m <- reshape2::melt(pseudo.d)

# Dot/Box plots
library(ggbeeswarm)
myCols <- readRDS("./data/annot/myColsCondition.rds")

#IFN
nnDotPlot(df = d.m,stimType = "IFN")
# LPS
nnDotPlot(df = d.m,stimType = "LPS")
# PMA
nnDotPlot(df = d.m,stimType = "PMA")

# Same, but for only Mono (Fig 4 in paper)

#IFN
IFNMono <- ggcustom(nnDotPlot(df = d.m,stimType = "IFN",cellType = "Mono"),clean = FALSE,splitLeg = FALSE)
# LPS
LPSMono <- ggcustom(nnDotPlot(df = d.m,stimType = "LPS",cellType = "Mono"),clean = FALSE,splitLeg = FALSE)
# PMA
PMAMono <- ggcustom(nnDotPlot(df = d.m,stimType = "PMA",cellType = "Mono"),clean = FALSE,splitLeg = FALSE)

cowplot::plot_grid(IFNMono,LPSMono,PMAMono,align="hv",nrow=1)


# Make density plot by Donor to show no donor-specific effect (Supp fig)
pseudo.d$Donor <- ATAC.meta$Donor
IFN.donor <- pseudo.d %>% filter(Condition %in% c("Control_1h","IFN_1h","IFN_6h","IFNGolgiPlug_6h") & cellType %in% "Mono") %>% ggplot(aes(x=IFN.scaled,color=Donor)) + geom_density() + facet_wrap(~Condition,nrow=1,scales = "free_y") + theme_classic() + scale_y_continuous(breaks = scales::pretty_breaks(n=2)) + theme(strip.background = element_rect(color=NA,fill=NA),legend.position = "none",strip.text = element_text(size=5)) + scale_color_manual(values=c("brown","skyblue","darkorange","darkorchid4"))
LPS.donor <- pseudo.d %>% filter(Condition %in% c("Control_1h","LPS_1h","LPS_6h","LPSGolgiPlug_6h")  & cellType %in% "Mono") %>% ggplot(aes(x=LPS.scaled,color=Donor)) + geom_density() + facet_wrap(~Condition,nrow=1,scales = "free_y") + theme_classic() + scale_y_continuous(breaks = scales::pretty_breaks(n=2)) + theme(strip.background = element_rect(color=NA,fill=NA),legend.position = "none",strip.text = element_text(size=5)) + scale_color_manual(values=c("brown","skyblue","darkorange","darkorchid4"))
PMA.donor <- pseudo.d %>% filter(Condition %in% c("Control_1h","PMA_1h","PMA_6h","PMAGolgiPlug_6h")  & cellType %in% "Mono") %>% ggplot(aes(x=PMA.scaled,color=Donor)) + geom_density() + facet_wrap(~Condition,nrow=1,scales = "free_y") + theme_classic() + scale_y_continuous(breaks = scales::pretty_breaks(n=2)) + theme(strip.background = element_rect(color=NA,fill=NA),legend.position = "none",strip.text = element_text(size=5)) + scale_color_manual(values=c("brown","skyblue","darkorange","darkorchid4"))
gComb <- cowplot::plot_grid(ggcustom(IFN.donor,splitLeg = FALSE,clean=FALSE),ggcustom(LPS.donor,splitLeg = FALSE,clean=FALSE),ggcustom(PMA.donor,splitLeg = FALSE,clean=FALSE),nrow=3)

gComb