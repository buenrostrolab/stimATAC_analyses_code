# Script for misc helper functions

### Author: Vinay Kartha
### Contact: <vinay_kartha@g.harvard.edu>
### Affiliation: Buenrostro Lab, Department of Stem Cell and Regenerative Biology, Harvard University

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


smoothScoresNN <- function(NNmat,
                           TSSmat,
                           geneList = NULL,
                           barcodesList=NULL,
                           nCores = 1)
{
  if (is.null(rownames(NNmat)))
    stop("NN matrix has to have matching cell IDs as rownames\n")
  if (!all.equal(rownames(NNmat), colnames(TSSmat)))
    stop("Nearest-neighbor matrix and TSS activity score matrix don't have matching cells ..\n")
  cat("Number of cells in supplied TSS matrix: ", ncol(TSSmat),
      "\n")
  cat("Number of genes in supplied TSS matrix: ", nrow(TSSmat),
      "\n")
  cat("Number of nearest neighbors being used per cell for smoothing: ",
      ncol(NNmat), "\n")
  if (!is.null(geneList)) {
    if (!(all(geneList %in% rownames(TSSmat)))) {
      cat("One or more of the gene names supplied is not present in the TSS matrix provided: \n")
      cat(geneList[!geneList %in% rownames(TSSmat)], sep = ", ")
      cat("\n")
      stop()
    }
    cat("Running TSS score smoothing for genes:", geneList,
        sep = "\n")
    cat("........\n")
    TSSmat <- TSSmat[rownames(TSSmat) %in% geneList, ]
  }
  else {
    if(nrow(TSSmat) > 10000){
      cat("Running smoothing for all genes in TSS matrix! (n = ",
          nrow(TSSmat), ") This is bound to take more time than querying specific markers ..\n",
          sep = "")
    }
  }
  opts <- list()
  pb <- txtProgressBar(min = 0, max = ncol(TSSmat), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  time_elapsed <- Sys.time()
  
  
  cl <- parallel::makeCluster(nCores)
  doSNOW::registerDoSNOW(cl)
  
  if(!is.null(barcodesList)){
    cat("Subsetting to ",length(barcodesList)," barcodes in dataset..\n")
    NNmat <- NNmat[barcodesList,]
  }
  cat("Running in parallel using ", nCores, "cores ..\n")
  matL <- foreach::foreach(x=1:nrow(NNmat),.options.snow = opts,.packages = c("Matrix","data.table","dplyr")) %dopar% {
    smoothedScore <- data.table(Matrix::rowMeans(TSSmat[, NNmat[x,]]))
    rownames(smoothedScore) <- rownames(TSSmat)
    colnames(smoothedScore) <- rownames(NNmat)[x]
    smoothedScore
  }
  
  parallel::stopCluster(cl)
  
  close(pb)
  cat("Merging results ..\n")
  smoothedMat <- dplyr::bind_cols(matL) %>% data.matrix() %>% Matrix(sparse=TRUE)
  rownames(smoothedMat) <- rownames(TSSmat)
  #stopifnot(all.equal(colnames(smoothedMat), colnames(TSSmat)))
  time_elapsed <- Sys.time() - time_elapsed
  cat(paste("\nTime Elapsed: ", time_elapsed, units(time_elapsed),
            "\n"))
  
  return(smoothedMat)
}

motifPeakZtest <- function(peakSet, ### peak_set: vector of peak ids you wish to test for motif enrichment
                           bgPeaks, ### bg_peaks: matrix of background peak selection iterations by chromVAR
                           tfMat ### tf: binary matrix of peak by motif
) {
  
  if(nrow(tfMat)!=nrow(bgPeaks))
    stop("Reference peak set used for TF and background peaks matrix must match..\n")
  
  if(!all(peakSet %in% 1:nrow(bgPeaks)))
    stop("One or more of the provided peak indices are out of the background peak set range ..\n")
  
  
  # get selected peak motif frequencies
  cat("Getting selected peak motif frequencies ..\n")
  
  # get frequency of motifs in test set (observed)
  p.tab <- Matrix::colMeans(tfMat[peakSet, ])
  
  # get the background frequency in peak sets of the same size
  cat("Getting background peak motif frequencies ..\n")
  # extract relevant rows (i.e. peakset being tested) from background peak matrix
  bg.f <- as.matrix(bgPeaks[peakSet, ])
  
  # calculate (background) motif frequencies in each iteration of background peaks corresponding to peakset
  bg.tab <- apply(bg.f[, c(1:ncol(bgPeaks))], 2, function(bg_iter) {
    
    b.i <- Matrix::colMeans(tfMat[bg_iter, ])
    return(b.i)
    
  })
  
  cat("Calculating empirical p values and z score p values ..\n")
  
  # loop over each motif and generate enrichment statistics compared to background
  m.p <- dplyr::bind_rows(lapply(names(p.tab), function(motif) {
    
    # calculate sd and mean frequencies for bg and selected peaks
    s <- sd(bg.tab[motif, ])
    bg_freq <- mean(bg.tab[motif, ])
    
    z_score <- (p.tab[motif] - bg_freq) / s
    
    # generate data.frame object of relevant statistics
    d <- data.frame(
      motifID = motif,
      gene = extractTFNames(motif),
      motif_obs_freq = p.tab[motif],
      motif_bg_freq = mean(bg.tab[motif, ]),
      motif_counts = p.tab[motif] * length(peakSet),
      emp_pval = 1 - (sum(bg.tab[motif, ] < p.tab[motif]) / ncol(bg.tab)),
      z_test = z_score,
      pval.z = 2 * pnorm(-abs(z_score)),
      signed.log10p = -log10(2 * pnorm(-abs(z_score))) * sign(z_score)
    )
    return(d)
  }))
  # sort by enrichment pval, motif observed frequency
  m.p <- dplyr::arrange(m.p,pval.z, motif_obs_freq)
  # return df of enrichment scores
  return(m.p)
}

extractTFNames <- function(motifIDs){
  if(all(grepl("_",motifIDs,fixed = TRUE))){
    sapply(strsplit(sapply(strsplit(motifIDs,"_LINE.",fixed=FALSE),"[[",2),"_",fixed=FALSE),"[[",2)
  } else {
    message("One or more provided motif IDs do not contain any '_' characters .. returning IDs as is")
    motifIDs
  }
}

all.unique <- function(x){
  length(x)==length(unique(x))
}


clean_theme <- function(){
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank())
}

ggcustom<- function(gg,clean=FALSE,splitLeg=TRUE,...){
  
  
  gg <- gg +  theme(plot.background=element_blank(),
                    panel.border = element_blank(),
                    axis.text = element_text(color="black",size=5.5),
                    axis.title = element_text(size=7),
                    line = element_line(size = 0.235), 
                    legend.title=element_text(size=7), 
                    legend.text=element_text(size=5),
                    ... # Additional optional parameters passed to theme()
  )
  
  if(clean)
    gg <- gg + clean_theme()
  
  if(splitLeg){
    leg <- cowplot::get_legend(gg)
    
    gg <- gg + theme(legend.position="none")
    
    print(gg)
    grid::grid.newpage()
    grid::grid.draw(leg) } else {
      print(gg)
    }
  
}

