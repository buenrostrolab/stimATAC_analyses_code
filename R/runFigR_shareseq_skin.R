# Script to run FigR on SHARE-seq skin data
# Data corresponding to https://www.sciencedirect.com/science/article/pii/S0092867420312538

### Author: Vinay Kartha
### Contact: <vinay_kartha@g.harvard.edu>
### Affiliation: Buenrostro Lab, Department of Stem Cell and Regenerative Biology, Harvard University


setwd("<data_analysis_folder>")


source("./code/FigR_functions.R")


setwd("/mnt/users/sai/skin.network")

atac.se <- readRDS("atac.se.rds")
SEsum_smooth <- readRDS("SEsum_smooth.rds")
RNA_nor <- readRDS("RNA_nor.rds") #unsmoothed
Cis.target.sig.filter <- readRDS("Cis.target.sig.filter.rds")
head(Cis.target.sig.filter)
cis.clean <- data.frame(Peak=Cis.target.sig.filter$index,
                        Gene=Cis.target.sig.filter$Gene,
                        rObs=Cis.target.sig.filter$Corr,
                        pvalZ=Cis.target.sig.filter$Pvalue)
cells <- colnames(atac.se)[atac.se$Type == "Trial60.skin"] 
DORC.k <- 20
figR.res <- runFigR(ATAC.se=atac.se[ ,cells], 
                    dorcTab = cis.clean,
                    genome="mm10",
                    dorcMat= as.matrix(SEsum_smooth[ ,cells]),
                    rnaMat=RNA_nor[ ,cells], # should smooth it...
                    dorcK=DORC.k, # How many dorc kNNs are we using to pool peaks
                    n_bg=50, # No. of background peaks to use for motif Z test
                    nCores = 10)

# All by all scatter

gAll <- ggplot(figR.res,aes(Corr.log10P,Enrichment.log10P,color=Score)) + 
  geom_point_rast(size=0.01,shape=16) + 
  theme_classic() + 
  scale_color_gradientn(colours = jdb_palette("solar_extra"),limits=c(-4,4),oob = scales::squish)
ggcustom(gAll,clean = FALSE,splitLeg = TRUE)

# Specific drivers

gScatter <- plotDrivers(figR.d = figR.res,marker = "Hoxc13")

ggcustom(gScatter,clean = FALSE,splitLeg = FALSE,legend.key.size = unit(0.3,"cm"))

# Heatmap
# Filtered set from Jason / Sai (see skin shiny app for txt file of cell type specific DORCs)
myDorcs <- readRDS("./data/TFnetwork/Sai_shareseq_skin_DORCS.rds")

share.heat <- plotfigRHeatmap(figR.d = figR.res,
                              score.cut = 2,
                              DORCs = myDorcs,
                              column_names_gp = gpar(fontsize=6),
                              show_row_dend = FALSE,
                              border=TRUE)

share.heat

# Mean rank plot
gBar <- rankDrivers(figR.d = figR.res,myLabels = c("Lef1","Dlx3","Grhl1","Gata6","Klf3","Hoxc13","Zeb1","Barx2","Pou2f3","Sox5"))
ggcustom(gBar,clean = FALSE,splitLeg = FALSE,legend.key.size = unit(0.3,"cm"))
