# Code to run FigR on Sai's share-seq skin data
# Data courtesy Sai

library(BuenRTools)
library(ggplot2)
library(reshape2)
library(BuenColors)
library(ggrepel)
library(circlize)

source("/mnt/Analyses/biorad/stim_06_17_2019/code/FigR_functions.R")
source("/mnt/Analyses/biorad/stim_06_17_2019/code/misc_helper_stim.R")


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
cells <- colnames(atac.se)[atac.se$Type == "Trial60.skin"] #[1:1000]
DORC.k <- 20
figR.res <- runFigR(ATAC.se=atac.se[ ,cells], 
                    dorcTab = cis.clean,
                    genome="mm10",
                    dorcMat= as.matrix(SEsum_smooth[ ,cells]),
                    rnaMat=RNA_nor[ ,cells], # should smooth it...
                    dorcK=DORC.k, # How many dorc kNNs are we using to pool peaks
                    n_bg=50, # No. of background peaks to use for motif Z test
                    nCores = 10)
setwd("/mnt/Analyses/biorad/stim_06_17_2019/")

saveRDS(figR.res, paste0("./data_freeze_10_16_2020/TFnetwork/SHAREseq_DORC_motif_enrich_k_",DORC.k,".rds"))

# All by all scatter

pdf("./data_freeze_10_16_2020/figures/DORCs/SHAREseq_FigR_all_scatter.pdf",height = 2,width = 2.3,useDingbats = FALSE)
gAll <- ggplot(figR.res,aes(Corr.log10P,Enrichment.log10P,color=Score)) + 
  geom_point_rast(size=0.01,shape=16) + 
  theme_classic() + 
  scale_color_gradientn(colours = jdb_palette("solar_extra"),limits=c(-4,4),oob = scales::squish)
ggcustom(gAll,clean = FALSE,splitLeg = TRUE)
dev.off()


# Specific drivers

gScatter <- plotDrivers(figR.d = figR.res,marker = "Hoxc13")

pdf(paste0("./data_freeze_10_16_2020/figures/figR/SHAREseq_FigR_Hoxc13_TF_regulation.pdf"),width=2.2,height=2.2,useDingbats = FALSE)
ggcustom(gScatter,clean = FALSE,splitLeg = FALSE,legend.key.size = unit(0.3,"cm"))
dev.off()

# Heatmap
# Filtered set from Jason / Sai (see skin shiny app for txt file of cell type specific DORCs, then filtered harder)
myDorcs <- readRDS("./data_freeze_10_16_2020/TFnetwork/Sai_shareseq_skin_DORCS.rds")

share.heat <- plotfigRHeatmap(figR.d = figR.res,
                              score.cut = 2,
                              DORCs = myDorcs,
                              column_names_gp = gpar(fontsize=6),
                              show_row_dend = FALSE,
                              border=TRUE)

pdf("./data_freeze_10_16_2020/figures/figR/SHAREseq_FigR_filt_heatmap.pdf",height = 5.5,width = 5)
share.heat
dev.off()

# Mean rank plot
gBar <- rankDrivers(figR.d = figR.res,myLabels = c("Lef1","Dlx3","Grhl1","Gata6","Klf3","Hoxc13","Zeb1","Barx2","Pou2f3","Sox5"))
pdf(paste0("./data_freeze_10_16_2020/figures/figR/SHAREseq_FigR_TF_DORC_regulation_all.pdf"),width=5,height=2.5,useDingbats = FALSE)
ggcustom(gBar,clean = FALSE,splitLeg = FALSE,legend.key.size = unit(0.3,"cm"))
dev.off()

# Network view
share.net <- plotfigRNetwork(figR.res,score.cut = 2,DORCs = myDorcs)
saveNetwork(share.net,file = "./data_freeze_10_16_2020/TFnetwork/SHAREseq_DORC_motif_enrich_k_20_filtNet.html",selfcontained = TRUE)
