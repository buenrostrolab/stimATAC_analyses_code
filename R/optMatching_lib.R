library(optmatch)
library(Matrix)
library(FNN)
library(dplyr)
  
# Function to plot resulting pairs
# Assumes condordant vectors for ATAC and RNA labels
# Also requires a data frame of UMAP coordinates with ATAC RNA labels as rownames
plotPairs <- function(ATAC,
                      RNA,
                      umap.df, # Assumes rownames are cell barcodes
                      max.show=300,
                      seed=123,
                      pairPoiintSize=1
                      ){
  
  cat("Assuming concordance in pairs ..\n")
  
  stopifnot(!is.null(rownames(umap.df)))
  stopifnot(all.unique(rownames(umap.df)))
  stopifnot(length(ATAC)==length(RNA))
  
  if(length(ATAC) > max.show){
    cat("Only highlighting",max.show,"pairs at random..\n")
    set.seed(seed)
    subsample <- sample(1:length(ATAC),max.show,replace = FALSE) 
    ATAC <- ATAC[subsample]
    RNA <- RNA[subsample]
  }
  
  d1 <- umap.df[ATAC,1:2]
  d2 <- umap.df[RNA,1:2]
  
  d1$pair <- factor(1:nrow(d1))
  d2$pair <- d1$pair
  
  d1$Group <- "ATAC"
  d2$Group <- "RNA"
  
  d <- rbind(d1,d2)
  
  
  require(ggrastr)
  ggplot(umap.df, aes(x = UMAP1, y = UMAP2)) +
    geom_point_rast(color="gray95",size=1.5)+
    geom_point(data=d,aes(UMAP1,UMAP2,color=Group),alpha=0.8,size=pairPoiintSize) + 
    scale_color_manual(values=c("cadetblue","darkorange"))+
    geom_line(data = d, aes(x = UMAP1, y = UMAP2, group = pair),alpha=0.4,color="darkgray")+
    theme_void()+
    labs(color="")+
    theme(legend.position = "bottom")+
    guides(colour = guide_legend(override.aes = list(size=3)))
}


# Function that evenly (randomly) samples larger dataset to smaller one
# Chunks up and returns list of chunked CCA ATAC/RNA components

chunk_CCA <- function(CCA_1,
                      CCA_2,
                      chunkSize=NULL,
                      seed=123){
  if(nrow(CCA_1) > nrow(CCA_2)){
    bigger <- CCA_1
    smaller <- CCA_2
  } else {
    bigger <- CCA_2
    smaller <- CCA_1
  }
  
  if(is.null(rownames(bigger)) | is.null(rownames(smaller)))
    stop("CCA matrices must have valid cell barcode names as rownames ..\n")
  
  # Shuffle them to make sure there's no bias in chunking
  set.seed(seed)
  bigger <- bigger[sample(1:nrow(bigger),nrow(bigger),replace = FALSE),]
  smaller <- smaller[sample(1:nrow(smaller),nrow(smaller),replace = FALSE),]
  
  cat("Number of cells in bigger dataset: ",nrow(bigger),"\n")
  cat("Number of cells in smaller dataset: ",nrow(smaller),"\n")
  cat("Difference in cell count between 2 datasets: ",nrow(bigger) - nrow(smaller),"\n")
  
  if(is.null(chunkSize))
     chunkSize <- nrow(smaller)
  
  cat("Chunking larger dataset to match smaller datset ..\n")
  cat("Chunk size n = ",chunkSize," cells\n")
  
  frac <- nrow(bigger) %/% chunkSize # Number of iterations of smaller in bigger
  remain <- nrow(bigger) %% chunkSize # Extra amount we have to sample after even chunking
  
  cat("Total number of chunks: ",length(1:frac) + sum(remain > 0),"\n")
  
  bigCells <- c()
  smallCells <- c()
  
  set.seed(seed)  
  chunkList <- list()
  
  for(i in 1:frac){
    # What has already not been sampled
    bigCellsBag <- rownames(bigger)[!rownames(bigger) %in% bigCells]
    smallCellsBag <- rownames(smaller) # Sample from all cells
    
    # Fetch chunk for each dataset
    bigCellsChunk <- sample(bigCellsBag,chunkSize,replace=FALSE)
    smallCellsChunk <- sample(smallCellsBag,chunkSize,replace=FALSE)
    
    # Update wha thas been sampled already
    bigCells <- c(bigCells,bigCellsChunk)
    smallCells <- c(smallCells,smallCellsChunk)
    
    # List of bigger and smaller CCA PCs for chunk
    # NOTE: Each chunk is always returned as Bigger 1st, smaller 2nd (ATAC, then RNA for us)
    chunkList[[i]] <- list(bigger[bigCellsChunk,],
                          smaller[smallCellsChunk,])
  }
  
  
  
  
  # Add extra cells for last uneven chunk (if any)
  if(remain > 0){
  chunkList[[i+1]] <- list(bigger[rownames(bigger)[!rownames(bigger) %in% bigCells],],
                           smaller[sample(rownames(smaller),remain,replace=FALSE),])
  }
  chunkList 
}

  
# Function that gives you an "up-sampled" dataset given two paired datasets
# Draws cells at random from smaller one till it matches larger one
upSample <- function(mat1,
                     mat2,
                     seed=123){
  if(nrow(mat1) > nrow(mat2)){
    bigger <- mat1
    smaller <- mat2 } else {
      bigger <- mat2
      smaller <- mat1
    }
  
  
  frac <- nrow(bigger) %/% nrow(smaller)
  rem <- nrow(bigger) %% nrow(smaller)
  
  cat("Difference in cell count between 2 datasets: ",nrow(bigger) - nrow(smaller),"\n")
  
  cat("Upsampling smaller to match larger ..\n")
  
  set.seed(seed)
  extraCells <- sample(1:nrow(smaller),rem,replace = FALSE)
  
  # Add chunk
  if(frac > 1)
    smaller <- smaller[rep(1:nrow(smaller),frac),]
  
  # Add extra
  smaller <- rbind(smaller,smaller[extraCells,])
  
  stopifnot(nrow(smaller)==nrow(bigger))
  
  cat("Both datasets now have ",nrow(smaller), " cells\n")
  
  if(length(intersect(rownames(smaller),rownames(mat1))) > 0) {
    list(smaller,bigger)  
  } else {
    list(bigger,smaller)  
  }
  
}


upSampleNew <- function(mat,
                     to,  
                     seed=123){

  
  to <- round(to)
  
  stopifnot(to > nrow(mat))
  
  frac <- to %/% nrow(mat)
  rem <- to %% nrow(mat)
  
  cat("Number of cells being sampled (added) : ",to - nrow(mat),"\n")
  
  cat("Upsampling smaller dataset ..\n")
  
  # Random extra cells
  set.seed(seed)
  extraCells <- sample(1:nrow(mat),rem,replace = FALSE)
  
  # Add chunk
  if(frac > 1)
    mat <- mat[rep(1:nrow(mat),frac),]
  
  # Add extra
  mat <- rbind(mat,mat[extraCells,])
  
  
  cat("Dataset now has ",nrow(mat), " cells\n")
  
  mat
}

# Function that gives you an "up-sampled" dataset given two paired datasets
# Draws cells using density-based upsampling from smaller one till it matches larger one
smartUpSample <- function(mat1,
                          mat2,
                          seed=123,
                          dist.use="pearson"){
  if(nrow(mat1) >= nrow(mat2)){
    bigger <- mat1
    smaller <- mat2 } else {
      bigger <- mat2
      smaller <- mat1
    }
  
  big.cells <- rownames(bigger)
  small.cells <- rownames(smaller)
  
  if(is.null(big.cells) | is.null(small.cells))
    stop("Both ATAC and RNA CCA matrices must have rownames (cell IDs) ..\n")
  
  # Perform Greedy pairing for these cells based on pearson correlation (max)
  cat("Performing initial greedy pairing using ",dist.use, "as metric ..\n")
  # Rows are bigger matrix cells, cols are smaller matrix cells
  cca.cor <- cor(t(bigger),t(smaller),method = dist.use)
  cor.max <- apply(cca.cor,1,function(x) names(which.max(x)))
  
  stopifnot(all.equal(rownames(cca.cor),names(cor.max)))
  
  pairing.greedy <- data.frame("bigger"=names(cor.max),"smaller"=as.character(cor.max),stringsAsFactors = FALSE)
  
  rm(cca.cor)
  gc()
  
  # E.g. For each RNA cell, see how many ATAC cells it paired to
  numBiggerpairs <- c()
  for(i in small.cells){
    numBiggerpairs <- c(numBiggerpairs,sum(pairing.greedy$smaller %in% i))
  }
  
  
  # Smooth the above number based on RNA CCA knn
  cat("Getting kNN of mapped neigbors ..\n")
  smaller.knn <- FNN::get.knn(smaller,k=200)$nn.index
  rownames(smaller.knn) <- small.cells
  
  numBiggerpairssmoothed <-  unlist(lapply(1:length(numBiggerpairs),function(i){
    mean(numBiggerpairs[smaller.knn[i,]])
  }))
  
  samplingDensity <- numBiggerpairssmoothed/length(big.cells)
  
  cat("Up-sampling smaller dataset ..\n")
  set.seed(seed)
  small.upCells <- sample(small.cells,size = length(big.cells),replace = TRUE,prob = samplingDensity)
  cat("Finished\n")
  
  
  smaller.up <- smaller[small.upCells,]
  
  stopifnot(nrow(smaller.up)==nrow(bigger))
  
  cat("Both datasets now have ",nrow(bigger), " cells\n")
  
  cat("% smaller dataset coverage up-on up-sampling: ",round((length(unique(small.upCells))/length(small.cells)) *100,2),"\n")
  
  if(length(intersect(rownames(smaller),rownames(mat1))) > 0) {
    list(smaller.up,bigger)  
  } else {
    list(bigger,smaller.up)  
  }
  
}


all.unique <- function(x){
  length(x)==length(unique(x))
}



# These two functions are only for Bio-Rad isolate data that had cell barcodes named a certain way
evaluatePairClass <- function(ATAC,RNA){
  freqs <- table(ATAC=getCellType(ATAC),RNA=getCellType(RNA)) 
  freqs / rowSums(freqs) 
}


getCellType <- function(cellNames){
  splitAndFetch(cellNames,"_",2)
}

plotCellFreq <- function(ATAC,RNA){
  cellType.orig.d <- rbind(as.data.frame(table(getCellType(ATAC))),as.data.frame(table(getCellType(RNA))))
  cellType.orig.d$Assay <- rep(c("ATAC","RNA"),each=5)
  
  ggplot(cellType.orig.d,aes(x=Var1,y=Freq,fill=Var1)) + 
    geom_bar(stat="identity",position="dodge") + 
    scale_fill_manual(values=myCols)+ 
    theme_classic()+
    facet_wrap(~Assay,nrow=1)+
    theme(strip.background = element_blank(),strip.text = element_text(size=12),axis.text = element_text(color="black"))+
    labs(fill="Experiment",y="Cell count",x="")
  
}

# Main matching function
# By default, will return a data frame of pairs
# Uses names of barcodes if present in input PC matrices (rownames)
# Otherwise uses relative indices of input PC matrices
optPair <- function(ATACpcs, # ATAC cells x PCs matrix (from CCA). Needs valid rownames
                    RNApcs, # RNA cells x PCs matrix (from CCA). Needs valid rownames
                    doPermutation=FALSE, # Test FDR on pair?
                    returnDistMat=FALSE, # Return the sparse distance matrix used for matching?
                    k=20, # k-NN parameter used for applying constraints on ATAC-RNA pairs
                    distMat=NULL, # If you already computed ATAC x RNA distance (any type of distance)
                    pairKNN=NULL, # If you've already computed ATAC x RNA, and RNA x ATAC KNNs (any way) in list
                    nCores=1,
                    sanityLoop=FALSE,
                    forceSymmetry=TRUE # ATAC and RNA have to be same dimensions
                    ){
  
  
  numATAC <- nrow(ATACpcs)
  numRNA <- nrow(RNApcs)
  
  if(forceSymmetry & numATAC!=numRNA){
    message("forceSymmetry set to TRUE")
    stop("Number of cells in both datasets needs to be the same ..\n")
  }
  
    if(is.null(rownames(ATACpcs)) | is.null(rownames(RNApcs))){
      message("Input matrices missing rownames (cell IDs) to use .. Using relative index instead \n")
      labelsExist <- FALSE # Use indices instead of cell IDs
    } else {
      labelsExist <- TRUE
      }
  
  time_elapsed <- Sys.time()
  
  if(is.null(pairKNN)){
  cat("Choosing min k required based on available neighbors and input parameter ..\n")
  k <- min(k,c(numATAC,numRNA))
    
    if(k > 300 | sanityLoop){
      options("optmatch_max_problem_size" = Inf)
      message("Warning: Very large k will lead to exceedingly large computation times .. \n")
    }
    
  cat("Getting knn graph between ATAC and RNA cells in PC space ..\n")
  cat("Restricting to ",k, "nearest neighbors ..\n")
  pairKNN12 <- FNN::get.knnx(data = RNApcs,query = ATACpcs,k =k)
  pairKNN21 <- FNN::get.knnx(data = ATACpcs,query = RNApcs,k =k)
  
  pairKNN <- list(ATAC=pairKNN12,RNA=pairKNN21)
  } else {
    cat("Using supplied KNN graph ..\n")
    stopifnot(c("ATAC","RNA"))
  }
  
  if(!is.null(distMat)){
    stopifnot(nrow(distMat)==numATAC & ncol(distMat)==numRNA)
    
    cat("Using supplied distance matrix ..\n")
    
    # Convert KNN to adjacency matrix
    cat("Getting KNN adjacency matrix ..\n")
    knn_adj_mtx <- array(0, dim = c(numATAC,numRNA))
    for (i in 1:numATAC) knn_adj_mtx[i,pairKNN$ATAC[i,]] <- 1 # Match. Won't work for FNN returned KNN (no $nn.index)
    for (i in 1:numRNA) knn_adj_mtx[pairKNN$RNA[i,],i] <- 1 # Match. Won't work for FNN returned KNN (no $nn.index) 
    
    # Scalar multiplication between matrices (will become 0 if not in union of KNNs)
    distMat <- knn_adj_mtx * distMat
    
    rownames(distMat) <- paste0("ATAC_",1:numATAC)
    colnames(distMat) <- paste0("RNA_",1:numRNA)
    
  } else {
  
  if(sanityLoop){
  # Make empty (sparse) matrix of ATAC (row) x RNA (col)
  distMat <- Matrix(0,nrow = numATAC,ncol=numRNA,sparse=TRUE) # We will replace the 0's as Inf later
  
  rownames(distMat) <- paste0("ATAC_",1:numATAC)
  colnames(distMat) <- paste0("RNA_",1:numRNA)
  
  # Update distances from knn graph only for nearest neighbors
  # ATAC NNs first
  cat("\nAssigning ATAC distances based on RNA NNs ..\n")
  for(i in 1:numATAC){
    distMat[i,pairKNN$ATAC$nn.index[i,]] <- pairKNN$ATAC$nn.dist[i,]
  }
  # RNA
  cat("Assigning RNA distances based on ATAC NNs ..\n")
  for(i in 1:numRNA){
    distMat[pairKNN$RNA$nn.index[i,],i] <- pairKNN$RNA$nn.dist[i,]
  }
  
   } else {
  
  # Parallelize
  stopifnot(nCores < 10)
  if(nCores > 1){  
  cat("Using ",nCores, " cores ..\n\n")   
  
  cat("\nAssigning distances based on ATAC - RNA NNs ..\n")
  
  tmpMat.l <- parallel::mclapply(1:numATAC,function(x) {
    vec <- rep(0,numRNA)
    vec[pairKNN$ATAC$nn.index[x,]] <- pairKNN$ATAC$nn.dist[x,]
    data.table::data.table(vec)
  },mc.cores = nCores)
  
  if(any(unlist(lapply(tmpMat.l,is.na))))
    stop("One or more distances returned NAs. Check input matrices\n")
  if(any(unlist(lapply(tmpMat.l,is.null))))
    stop("One or more distances returned NULL. Check input matrices or for parallel failure\n")
  
  tmpMat <- Matrix::t(dplyr::bind_cols(tmpMat.l) %>% data.matrix() %>% 
    Matrix::Matrix(sparse=TRUE))
  
  cat("Assigning distances based on RNA - ATAC NNs ..\n")
  
  distMat.l <- parallel::mclapply(1:numRNA,function(x,M=tmpMat) {
    vec <- as.numeric(M[,x])
    vec[pairKNN$RNA$nn.index[x,]] <- pairKNN$RNA$nn.dist[x,]
    data.table::data.table(vec)
  },mc.cores = nCores)

  if(any(unlist(lapply(distMat.l,is.na))))
    stop("One or more distances returned NAs. Check input matrices\n")
  if(any(unlist(lapply(distMat.l,is.null))))
    stop("One or more distances returned NULL. Check input matrices or for parallel failure\n")
  
  cat("Merging parallel chunks ..\n")
  distMat <- dplyr::bind_cols(distMat.l) %>% data.matrix() %>% 
                  Matrix::Matrix(sparse = TRUE)
  
  } else {
    # Use regular lapply
    cat("\nAssigning distances based on ATAC - RNA NNs ..\n")
    
    tmpMat.l <- lapply(1:numATAC,function(x) {
      vec <- rep(0,numRNA)
      vec[pairKNN$ATAC$nn.index[x,]] <- pairKNN$ATAC$nn.dist[x,]
      data.table::data.table(vec)
    })
    
    if(sum(unlist(lapply(tmpMat.l,is.na)))!=0)
      stop("One or more distances returned NAs. Check input matrices\n")
    
    tmpMat <- Matrix::t(dplyr::bind_cols(tmpMat.l) %>% data.matrix() %>% 
                          Matrix::Matrix(sparse=TRUE))
    
    cat("Assigning distances based on RNA - ATAC NNs ..\n")
    
    distMat.l <- lapply(1:numRNA,function(x,M=tmpMat) {
      vec <- as.numeric(M[,x])
      vec[pairKNN$RNA$nn.index[x,]] <- pairKNN$RNA$nn.dist[x,]
      data.table::data.table(vec)
    })
    
    if(sum(unlist(lapply(distMat.l,is.na)))!=0)
      stop("One or more distances returned NAs. Check input or KNN matrices\n")
    
    distMat <- dplyr::bind_cols(distMat.l) %>% data.matrix() %>% 
      Matrix::Matrix(sparse = TRUE)
    
    } # End nCores 1 normal lapply
      # For apply cases, assign row/colnames after
     rownames(distMat) <- paste0("ATAC_",1:numATAC)
     colnames(distMat) <- paste0("RNA_",1:numRNA) 
     
   } # End not sanity loop case
  } # End if supplied own distance matrix
  
  # Make 0 entries infinite (distance of invalid pairs)
  distMat[distMat==0] <- Inf
  
  cat("\nDeterming pairs through optimized bipartite matching ..\n")
  myMatches <- suppressWarnings(optmatch::pairmatch(optmatch::as.InfinitySparseMatrix(as.matrix(distMat))))
  myMatches <- sort(myMatches) # Sort to get ATAC, RNA tuples
  
  
  if(length(myMatches)==0)
    stop("Matches could not be found .. Perhaps try adjusting the constraints to allow optimal matching to be solved?\n")
  
  if(any(is.na(myMatches)))
    warning("NA pairs exist ..\n")
  
  
  # Make sure pair groupings (factors) are adjacent
  stopifnot(all.equal(myMatches[seq(1,length(myMatches),2)],myMatches[seq(2,length(myMatches),2)],check.attributes=FALSE))
  
  
  # Sometimes this is arranged as RNA/ATAC instead of ATAC first
  # Check which is first and fetch corresponding labels
  ATACtups <- which(splitAndFetch(names(myMatches),"_",1) %in% "ATAC")[1]
  RNAtups <- ifelse(ATACtups==1,2,1)
  
  ATACp <- names(myMatches)[seq(ATACtups,length(myMatches),2)] # Names of ATAC cells
  RNAp <- names(myMatches[seq(RNAtups,length(myMatches),2)]) # Names of RNA cells
  
  # Checking if ATAC and RNA tuples are concordant
  stopifnot(all(splitAndFetch(ATACp,"_",1) %in% "ATAC"))
  stopifnot(all(splitAndFetch(RNAp,"_",1) %in% "RNA"))
  
  # This is just to make sure 1-1, with the names we gave 
  # (can still comprise actual doublets from upsampling if any)
  stopifnot(all.unique(ATACp) & all.unique(RNAp)) 
  
  # Get corresponding index relative to input matrix order
  ATACi <- as.numeric(splitAndFetch(ATACp,"_",2))
  RNAi <- as.numeric(splitAndFetch(RNAp,"_",2))
  
  myMatches.m <- matrix(c(ATACi,RNAi),ncol=2,byrow = FALSE) # ATAC col1, RNA col2
  
  cost <- distMat[myMatches.m] # Fetches distance per ATAC/RNA pair (in corresponding pair order)
  

  cat("Assembling pair list ..\n")
  # Make data frame of matches
  
  pairDF <- data.frame("ATAC"=ATACi,
                       "RNA"=RNAi,
                       "Distance"=cost)
  
  cat("Finished!\n")
  time_elapsed <- Sys.time() - time_elapsed
  
  cat(paste("Total Run-time: ",round(time_elapsed,2),units(time_elapsed),"\n\n"))
  
  permuted.distances <- c()
  if(doPermutation){
    cat("Running permutations ..\n")
    numPerm <- 100
    for(i in 1:numPerm){
      set.seed(i) 
      cat("iteration ",i," ..\n")
      permuted.distances <- cbind(permuted.distances,distMat[matrix(c(ATACi,sample(unique(c(pairKNN$ATAC$nn.index)),length(RNAi),replace = FALSE)),ncol=2,byrow = FALSE)]) # ATAC col1, RNA col2
    }
    
    perm_p <- rowSums(permuted.distances < cost) / numPerm
    pairDF$perm_pval <- perm_p
    
  }
  
  pairDF <- pairDF %>% arrange(ATAC)
    
  # Convert to labels if they exist
  if(labelsExist){
    pairDF$ATAC <- rownames(ATACpcs)[pairDF$ATAC] # ATAC cell labels
    pairDF$RNA <- rownames(RNApcs)[pairDF$RNA] # RNA cell labels
  }
  
  # Assemble resulting list
  pairList <- list()
  pairList$pairs <- pairDF # Pair data frame
  pairList$Kparam <- k # K NN param used
  
  if(!labelsExist){
  cat("% of RNA cell barcodes accounted for in match : ",round((sum(1:nrow(RNApcs) %in% pairList$pairs$RNA)/nrow(RNApcs)*100),2), "% \n")
  cat("% of ATAC cell barcodes accounted for in match : ", round((sum(1:nrow(ATACpcs) %in% pairList$pairs$ATAC)/nrow(ATACpcs)*100),2), "% \n\n")
  } else {
  cat("% of RNA cell barcodes accounted for in match : ",round((sum(rownames(RNApcs) %in% pairList$pairs$RNA)/nrow(RNApcs)*100),2), "% \n")
  cat("% of ATAC cell barcodes accounted for in match : ", round((sum(rownames(ATACpcs) %in% pairList$pairs$ATAC)/nrow(ATACpcs)*100),2), "% \n\n")
  }
  
  cat("Frequencies of RNA barcode supermatches : \n")
  print(table(table(pairList$pairs$RNA)))
  
  # Return the pairing and the sparse-ified distance matrix used as input to pairing, and the k param
  if(returnDistMat)
    pairList$dist <- distMat # Note this won't have the labels as row/column names since they might not be unique (just has "ATAC_1" etc.)
  
  pairList
  
}
