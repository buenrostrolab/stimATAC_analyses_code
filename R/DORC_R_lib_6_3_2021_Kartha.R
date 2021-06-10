splitAndFetch <- function(vec,
                          delim,
                          part){
  if(length(part)==1){
    sapply(strsplit(as.character(vec),delim,fixed=TRUE),"[[",part) } else {
      sapply(strsplit(as.character(vec),delim,fixed = TRUE),function(x) paste(x[part],collapse = delim))
    }
}

centerCounts <- function(obj,
                         doInChunks=TRUE,
                         chunkSize=1000){
  if(!class(obj) %in% c("SummarizedExperiment","RangedSummarizedExperiment","dgCMatrix","dgeMatrix","Matrix"))
    stop("Supplied object must be either of class SummarizedExperiment or sparse Matrix ..\n")
  
  if(ncol(obj) > 10000)
    doInChunks <- TRUE
  
  if(doInChunks){
    cat("Centering counts for cells sequentially in groups of size ",
        chunkSize, " ..\n\n")
    starts <- seq(1,ncol(obj),chunkSize)
  } else{
    starts <- 1
  }
  
  counts.l <- list()
  
  for(i in 1:length(starts)){
    beginning <- starts[i]
    if(i==length(starts)) # If it's the last one
    {
      ending <- ncol(obj)
    } else {
      ending <- starts[i]+chunkSize-1
    }
    
    cat("Computing centered counts for cells: ",beginning," to ", ending,"..\n")
    
    if(class(obj) == "RangedSummarizedExperiment" | class(obj)=="SummarizedExperiment"){
      m <- SummarizedExperiment::assay(obj[, beginning:ending])} else {
        m <- obj[,beginning:ending] # Assumes Matrix format
      }
    cellMeans <- Matrix::colMeans(m)
    cat("Computing centered counts per cell using mean reads in features ..\n\n")
    # Center cell counts based on its mean RIP count
    cCounts <- Matrix::t(Matrix::t(m)/cellMeans)
    
    counts.l[[i]] <- cCounts
    
    gc()
  }
  
  cat("Merging results..\n")
  centered.counts <- do.call("cbind",counts.l)
  cat("Done!\n")
  
  if(class(obj) == "RangedSummarizedExperiment" | class(obj)=="SummarizedExperiment"){
    SummarizedExperiment::assay(obj) <- centered.counts
    return(obj)
  } else {
    return(centered.counts)
  }
}

chunkCore <- function(chunk,
                      A, # ATAC matrix
                      R, # RNA matrix
                      O, # Gene-Peak overlap pairing data.frame
                      met # Correlation method ("spearman" or "pearson")
){
  
  # Get indices of genes and peaks from overlap object for chunk
  # Assumes query hits are genes and subject hits are peaks in the overlap object
  geneIndices <- O$Gene[chunk[1]:chunk[2]]
  peakIndices <- O$Peak[chunk[1]:chunk[2]]
  
  pairnames <- cbind(rownames(A)[peakIndices],rownames(R)[geneIndices])
  
  uniquegenes <- unique(geneIndices)
  uniquepeaks <- unique(peakIndices)
  
  
  M1 <- as.matrix(Matrix::t(A[uniquepeaks,,drop=FALSE])) # In case only 1 match, keep matrix structure
  M2 <- as.matrix(Matrix::t(R[uniquegenes,,drop=FALSE])) # In case only 1 match, keep matrix structure
  
  # Peak x Gene correlation matrix, subset by peak-gene pair names to get corresponding correlation vector
  # NOTE: This correlation call fails if you have maps with just 1 gene / peak. This is unlikely for large chunk sizes
  cor(x = M1,y = M2,method = met)[pairnames]
  
}

PeakGeneCor <- function(ATAC, # Normalized reads in peaks counts (rownames should  be "Peak1","Peak2" etc.)
                        RNA, # Normalized gene expression counts
                        OV, # Gene TSS - Peak overlap pairs object (Genes: query, Peaks: subject)
                        ncores=4,
                        chunkSize=1000,
                        metric="spearman",
                        bg=NULL){
  
  stopifnot(ncol(ATAC)==ncol(RNA))
  
  if(chunkSize > 1000)
    stop("Do not specify very large chunk sizes. Please use chunkSize < 1000")
  
  
  # Number of total gene-peak pairs to chunk up for parallelization
  n <- length(OV)
  starts <- seq(1, n, chunkSize)
  ends <- starts + chunkSize -1
  ends[length(ends)] <- n
  
  OVd <- OV %>% as.data.frame() %>% rename("Gene"="queryHits","Peak"="subjectHits")
  
  chunkList <- mapply(c, starts, ends, SIMPLIFY = FALSE)
  
  time_elapsed <- Sys.time()
  
  cat("Running in parallel using ", ncores, "cores ..\n")
  
  cat("Computing observed correlations ..\n")
  
  corList <- pbmcapply::pbmclapply(X=chunkList,
                                   FUN=function(x) {chunkCore(chunk=x,A=ATAC,R=RNA,O=OVd,met=metric)},mc.cores = ncores)
  
  
  if(any(unlist(sapply(corList,is.null)))){
    message("One or more of the chunk processes failed unexpectedly (returned NULL) ..")
    message("Please check to see you have enough cores/memory allocated")
    message("Also make sure you have filtered down to non-zero peaks/genes")
  }
  
  OVd$rObs <- unlist(corList)
  
  
  cat("Finished!\n")
  
  time_elapsed <- Sys.time() - time_elapsed
  cat(paste("\nTime Elapsed: ",time_elapsed, units(time_elapsed)),"\n\n")
  
  if(!is.null(bg)){
    n_iter <- ncol(bg)
    cat("Computing background correlations ..\n")
    
    time_elapsed <- Sys.time()
    
    bgCor <- foreach(i=1:n_iter,.combine = 'cbind',
                     .export = c("chunkCore","t"),.packages = c("pbmcapply","Matrix")) %do% {
                       OVdBg <- OVd[,1:2] # Initialize gene-peak pairing to observed
                       OVdBg$Peak <- bg[OVdBg$Peak,i] # Swap actual peaks with bg peaks for given iteration in pairing
                       bgCorList <- pbmcapply::pbmclapply(X=chunkList,
                                                          FUN=function(x) {chunkCore(chunk=x,A=ATAC,R=RNA,O=OVdBg,met=metric)},mc.cores = ncores)
                       unlist(bgCorList) # Vector of background permuted correlation values for that set of background peaks
                     }
    
    
    if(sum(is.null(bgCor))!=0 | sum(is.na(bgCor))!=0)
      stop("One or more of the chunk processes failed unexpectedly (returned NULL) .. Please check to see you have enough cores/m
           emory allocated")
    
    time_elapsed <- Sys.time() - time_elapsed
    cat(paste("\nTime Elapsed: ",time_elapsed, units(time_elapsed)),"\n\n")
    
    colnames(bgCor) <- paste0("rBg",1:ncol(bgCor))
    OVd <- cbind(OVd,bgCor)
  }
  
  return(OVd)
}

#' Gene-Peak correlations used for calling domains of regulatory chromatin (DORCs)
#'
#'Function to compute correlation between RNA expression and peak accessibility for peaks falling within a window around each gene, across single cells, using background peak correlations for significance testing
#'@param ATAC.se SummarizedExperiment object of the scATAC-seq reads in peak counts.
#'@param RNAmat Matrix object of the normalized scRNA-seq counts. For example, this will be the count data housed within the `@assays$RNA@data` field if processed using Seurat. Must have same number of cells (i.e. matched) with the scATAC-seq data
#'@param genome character specifying the reference genome build to use. Must be one of "hg19", "hg38" or "mm10", with no default
#'@param TSSwindow numeric specifying the window size (in base pairs) to pad around either side of each TSS, to fetch overlapping peaks. Default is 50 kb (so 100 kb window is drawn arund each TSS)
#'@param normalizeATACmat boolean indicating whether or not to normalize the counts present in the ATAC.se object prior to computing correlations. Default is TRUE (i.e. assumes peak counts in ATAC.se are raw).
#'@param n_bg numeric indicating the number of background correlations to compute per gene-peak pair. Default is 100 (increasing this can significantly increase run time)
#'@param p.cut numeric indicating p-value cut-off to apply to gene-peak correlations. Default is NULL (i.e. don't apply any filtering of results). Good option is to set to 0.05.
#'@param nCores numeric indicating the number of cores to use if parallelizing tasks
#'@return a data.frame of correlations for each gene-peak overlap pair, with the observed correlation level and significance p-value based on background correlations per pair
#'@import Matrix dplyr pbmcapply SummarizedExperiment GenomicRanges
#'@export
#'@author Vinay Kartha
runGenePeakcorr <- function(ATAC.se,
                            RNAmat,
                            genome, # Must be one of "hg19", "mm10", or "hg38"
                            geneList=NULL, # 2 or more valid gene symbols
                            windowPadSize=50000,
                            normalizeATACmat=TRUE,
                            nCores=4,
                            n_bg=100,
                            p.cut=NULL # Optional, if specified, will only return sig hits
) {
  
  stopifnot(inherits(ATAC.se,"RangedSummarizedExperiment"))
  stopifnot(inherits(RNAmat,c("Matrix","matrix")))
  
  if(!all.equal(ncol(ATAC.se),ncol(RNAmat)))
    stop("Input ATAC and RNA objects must have same number of cells")
  
  message("Assuming paired scATAC/scRNA-seq data ..")
  
  # Function needs rownames for both matrices or gives error
  rownames(ATAC.se) <- paste0("Peak",1:nrow(ATAC.se))
  ATACmat <- assay(ATAC.se) # Rownames preserved
  
  # Normalize peak counts
  if(normalizeATACmat)
    ATACmat <- centerCounts(ATACmat) # Rownames preserved
  
  if(is.null(rownames(RNAmat)))
    stop("RNA matrix must have gene names as rownames")
  
  # Check for peaks/genes with 0 accessibility/expression
  
  if(any(Matrix::rowSums(assay(ATAC.se))==0)){
    message("Peaks with 0 accessibility across cells exist ..")
    message("Removing these peaks prior to running correlations ..")
    peaksToKeep <- Matrix::rowSums(assay(ATAC.se))!=0
    ATAC.se <- ATAC.se[peaksToKeep,] # Subset ranges
    ATACmat <- ATACmat[peaksToKeep,]
    message("Important: peak indices in returned gene-peak maps are relative to original input SE")
  }
  
  
  peakRanges <- granges(ATAC.se) # Peak ranges
  
  if(any(Matrix::rowSums(RNAmat)==0)){
    message("Genes with 0 expression across cells exist ..")
    message("Removing these genes prior to running correlations ..")
    genesToKeep <- Matrix::rowSums(RNAmat)!=0
    RNAmat <- RNAmat[genesToKeep,]
  }
  
  cat("Number of peaks in ATAC data:",nrow(ATACmat),"\n")
  cat("Number of genes in RNA data:",nrow(RNAmat),"\n")
  
  
  if (!genome %in% c("hg19", "hg38", "mm10"))
    stop("You must specify one of hg19, hg38 or mm10 as a genome build for currently supported TSS annotations..\n")
  switch(genome, hg19 = {
    TSSg <- hg19TSSRanges
  }, hg38 = {
    TSSg <- hg38TSSRanges
  }, mm10 = {
    TSSg <- mm10TSSRanges
  })
  
  # Keep genes that have annotation and are in RNA matrix
  names(TSSg) <- as.character(TSSg$gene_name)
  
  if(!is.null(geneList)){
    if(length(geneList)==1)
      stop("Please specify more than 1 valid gene symbol")
    
    if(any(!geneList %in% names(TSSg))){
      cat("One or more of the gene names supplied is not present in the TSS annotation specified: \n")
      cat(geneList[!geneList %in% names(TSSg)], sep = ", ")
      cat("\n")
      stop()
    }
    
    TSSg <- TSSg[geneList]
  }
  
  # Checking in case some genes in RNA don't overlap our TSS annotations
  genesToKeep <- intersect(names(TSSg),rownames(RNAmat))
  cat("\nNum genes overlapping TSS annotation and RNA matrix being considered: ",length(genesToKeep),"\n")
  
  # Match gene order in RNA matrix and TSS ranges
  RNAmat <- RNAmat[genesToKeep,]
  TSSg <- TSSg[genesToKeep]
  
  # Pad TSS by this much *either side*
  TSSflank <- GenomicRanges::flank(TSSg,
                                   width = windowPadSize,
                                   both = TRUE)
  
  # Get peak summit
  cat("\nTaking peak summits from peak windows ..\n")
  peakSummits <- resize(peakRanges,width = 1,fix = "center")
  
  # Find overlap of all peaks to all genes given window
  # Subject is Peaks, query is Gene
  cat("Finding overlapping peak-gene pairs ..\n")
  genePeakOv <- findOverlaps(query = TSSflank,subject = peakSummits)
  numPairs <- length(genePeakOv)
  
  cat("Found ",numPairs,"total gene-peak pairs for given TSS window ..\n")
  
  cat("Number of peak summits that overlap any gene TSS window: ",length(unique(subjectHits(genePeakOv))),"\n")
  cat("Number of gene TSS windows that overlap any peak summit: ",length(unique(queryHits(genePeakOv))),"\n\n")
  
  # For each gene, determine observed correlation of each overlapping peak to its associated gene (gene expression)
  
  # For each of those genes, also determine correlation based on background peaks (run in parallel) and save
  # Do this in parallel, and get data frame of gene-peak-pearson values
  # Fetch background peaks for each peak tested (i.e. that has overlap in window with gene)
  set.seed(123)
  cat("Determining background peaks ..\n")
  
  if(is.null(rowData(ATAC.se)$bias)){
    if(genome %in% "hg19")
      myGenome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
    if(genome %in% "mm10")
      myGenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
    if(genome %in% "hg38")
      myGenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    
    ATAC.se <- chromVAR::addGCBias(ATAC.se,genome=myGenome) }
  
  cat("Using ",n_bg," iterations ..\n\n")
  
  set.seed(123)
  bg <- chromVAR::getBackgroundPeaks(ATAC.se,niterations=n_bg)
  
  cat("Computing gene-peak correlations ..\n")
  
  pairsPerChunk <- 1000
  
  # This defines the outer (larger chunks)
  largeChunkSize <- 20000
  
  startingPoint <- 1 # If for any reason the workers fail, resume from where it failed by specifying the starting point here
  chunkStarts <- seq(startingPoint, numPairs, largeChunkSize)
  chunkEnds <- chunkStarts + largeChunkSize -1
  chunkEnds[length(chunkEnds)] <- numPairs
  
  library(doParallel)
  
  dorcList <- list()
  for(i in 1:length(chunkStarts)){
    cat("Running pairs: ",chunkStarts[i], "to",chunkEnds[i],"\n")
    # This fill further chunk this up and run in parallel, saving the merged output ObsCor
    ObsCor <- PeakGeneCor(ATAC = ATACmat,
                          RNA = RNAmat,
                          OV = genePeakOv[chunkStarts[i]:chunkEnds[i]],
                          chunkSize = pairsPerChunk,
                          ncores = nCores,
                          bg = bg)
    gc()
    
    dorcList[[i]] <- ObsCor
  }
  
  cat("\nMerging results ..\n")
  dorcTab <- bind_rows(dorcList)
  
  cat("Performing Z-test for correlation significance ..\n")
  permCols <- 4:(ncol(bg)+3)
  
  
  # Filter to positive correlations
  cat("Only considering positive correlations ..\n")
  dorcTabFilt <- dorcTab %>% filter(rObs > 0)
  
  # Remove multi-mapping peaks (force 1-1 mapping)
  cat("Keeping max correlation for multi-mapping peaks ..\n")
  dorcTabFilt <- dorcTabFilt %>% group_by(Peak) %>% filter(rObs==max(rObs))
  
  
  # Swap gene number for gene symbol from TSS annotation lookup
  dorcTabFilt$Gene <- as.character(TSSg$gene_name)[dorcTabFilt$Gene]
  
  # Swap peak numbers to match reference input peak numbers
  # This only changes if some peaks had zero accessibility and were filtered out internally
  # Use rownames from reference matching
  dorcTabFilt$Peak <- as.numeric(splitAndFetch(rownames(ATACmat)[dorcTabFilt$Peak],"Peak",2))
  
  # # Z test pval
  dorcTabFilt$rBgSD <- matrixStats::rowSds(as.matrix(dorcTabFilt[,permCols]))
  dorcTabFilt$rBgMean <- rowMeans(dorcTabFilt[,permCols])
  dorcTabFilt$pvalZ <- 1-stats::pnorm(q = dorcTabFilt$rObs,mean = dorcTabFilt$rBgMean,sd = dorcTabFilt$rBgSD)
  
  
  cat("\nFinished!\n")
  
  if(!is.null(p.cut)){
    cat("Using significance cut-off of ",p.cut," to subset to resulting associations\n")
    dorcTabFilt <- dorcTabFilt[dorcTabFilt$pvalZ <= p.cut,] # Subset to significant correlations only
  }
  
  return(as.data.frame(dorcTabFilt[,c("Peak","Gene","rObs","pvalZ")],stringsAsFactors=FALSE))
}