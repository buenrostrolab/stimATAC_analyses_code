# Script to run core FigR function on stimulation data using new human motif database

### Author: Vinay Kartha
### Contact: <vinay_kartha@g.harvard.edu>
### Affiliation: Buenrostro Lab, Department of Stem Cell and Regenerative Biology, Harvard University

setwd("<data_analysis_folder>")


source("./code/FigR_functions.R")

# ---------------------------------------------------------------------------------- Load stim data / DORC calls
# Load SE object
SE.filt <- readRDS("./data/SE/ATAC_stim_paired.rds")
SE.filt

# DORC calls, all cells
sigGP <- read.table("./data/DORCs/DORC_calls_stim_allCells_0.05.tsv",header=TRUE,stringsAsFactors = FALSE,sep="\t")

# # Only work in space of DORC genes
dorcGenes <- readRDS("./data/DORCs/dorcGenes.rds")
sigGP <- sigGP[sigGP$Gene %in% dorcGenes,]

dim(sigGP)

# Normalized counts from SEs
# Load paired RNA
RNA.SE <- readRDS("./data/SE/RNA_stim_paired.rds")
names(assays(RNA.SE))
rnaMat <- assays(RNA.SE)$counts.norm
rownames(rnaMat) <- gsub(x=rownames(rnaMat),pattern = "-",replacement = "",fixed = TRUE)

DORC.SE <- readRDS("./data/SE/ATAC_stim_DORC_prom_enhancer_scores.rds")
names(assays(DORC.SE))

DORC.SE <- DORC.SE[,colnames(SE.filt)] # Subset to paired cells only

dorcMat <- assays(DORC.SE)$P + assays(DORC.SE)$E 
dim(dorcMat)

# ---------------------------------------------------------------------------------- Smooth data
lsi <- read.table("./data/ArchR/scATAC_stim_LSI.tsv",sep="\t",header=TRUE,stringsAsFactors=FALSE,row.names = 1)
head(lsi)

stopifnot(all(colnames(DORC.SE) %in% rownames(lsi)))

# Subset to paired ATAC
lsi <- lsi[colnames(DORC.SE),]

# Get KNNs
lsi.knn <- FNN::get.knn(lsi,k=30)$nn.index
rownames(lsi.knn) <- colnames(DORC.SE)

library(doParallel)
dorcMat.smoothed <- smoothScoresNN(lsi.knn,dorcMat,nCores = 6)
stopifnot(all.equal(rownames(dorcMat.smoothed),dorcGenes))

rm(DORC.SE)
gc()

# Just so that the smoothing function will work, since it checks for matching attributes
rownames(lsi.knn) <- colnames(RNA.SE)

# Run only on TFs to save time
human_pwms_v3 <- readRDS("./data/cisBP_human_pfms_2021.rds")
rnaMat.smoothed <- smoothScoresNN(NNmat = lsi.knn,TSSmat = rnaMat,nCores = 6,geneList=intersect(rownames(rnaMat),names(human_pwms_v3)))
gc()

# ---------------------------------------------------------------------------------- Run FigR

stim_FigR <- runFigR(ATAC.se = SE.filt,
                           dorcK = 30,
                           dorcTab = sigGP,
                           genome = "hg19", 
                           dorcMat = dorcMat.smoothed,
                           rnaMat = rnaMat.smoothed,
                           n_bg = 50,
                           nCores = 8)

library(ggplot2)
library(ggrastr)
library(BuenColors)
gAll <- ggplot(stim_FigR,aes(Corr.log10P,Enrichment.log10P,color=Score)) + 
  geom_point_rast(size=0.01,shape=16) + 
  theme_classic() + 
  scale_color_gradientn(colours = jdb_palette("solar_extra"),limits=c(-4,4),oob = scales::squish)

# Mean plot
library(ggrepel)
rankDrivers(figR.d = stim_FigR) 

# Scatter plot of drivers associated with a given DORC
plotDrivers(figR.d = stim_FigR,marker = "MX1")
plotDrivers(figR.d = stim_FigR,marker = "REL")
plotDrivers(figR.d = stim_FigR,marker = "TRAF1")

# Load DORC-implicated GWAS SNP P mat
SNPPmat <- readRDS("./data/GWAS/GWAS_DORC_SNPoV_unfiltered_P_mat.rds")
SNPDORCs <- rownames(SNPPmat)  # DORCs of interest
length(SNPDORCs)

# Heatmap 
figR_heat <- plotfigRHeatmap(figR.d = stim_FigR,
                             score.cut = 1.5,
                             DORCs = SNPDORCs,
                             column_names_gp=gpar(fontsize=5),
                             show_row_dend = FALSE,
                             row_names_side = "left")

# D3Network
plotfigRNetwork(stim_FigR,score.cut = 1.5,DORCs = SNPDORCs)