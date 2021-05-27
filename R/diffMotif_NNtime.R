library(SummarizedExperiment)
library(chromVAR)
library(BuenRTools)
library(ggrepel)
library(BuenColors)


# Script to help identify what the early responders are, and score cells by nearest neighbor time points

setwd("/mnt/Analyses/biorad/stim_06_17_2019/")

# source("./code/myggTheme.R")
source("./code/misc_helper_stim.R")

# Load ATAC SE (paired)
# NOTE that we load the paired dataset here
ATAC.SE <- readRDS("./data_freeze_10_16_2020/SE/ATAC_stim_paired.rds")
ATAC.meta <- as.data.frame(colData(ATAC.SE),stringsAsFactors=FALSE)


# Script to do differential testing on leaders and laggers (and other groupings) using nn stim time

# V2 below is for weights 0-2; 
# See computeNNtime.R for how this was generated
pseudoFile <- "./data_freeze_10_16_2020/nnTime/LSI_NN_stim_pseudotimev2.tsv"

pseudo.d <- read.table(pseudoFile,sep="\t",header=TRUE,stringsAsFactors = FALSE) 

identical(rownames(ATAC.meta),rownames(pseudo.d))

# Focus on LPS 1h/6h Mono cells
stim <- "LPS"
cellType <- "Mono"

# 1h and 6h Mono cells for the stim in question
ctrl1hMono <- ATAC.meta$pairedLabel2 %in% cellType & ATAC.meta$Condition %in% "Control_1h"
stim1hMono <- ATAC.meta$pairedLabel2 %in% cellType & ATAC.meta$Condition %in% paste(stim,"1h",sep="_") 
stim6hMono <- ATAC.meta$pairedLabel2 %in% cellType & ATAC.meta$Condition %in% paste(stim,"6h",sep="_")

table(ctrl1hMono)
table(stim1hMono)
table(stim6hMono)

# For footprinting and other analyses (if using)
write.table(rownames(ATAC.meta)[stim1hMono],paste0("./footprint/1h_vs_6h/",stim,"/",stim,"_1h_Mono.tsv"),sep="\n",row.names=FALSE,col.names = FALSE,quote=FALSE)
write.table(rownames(ATAC.meta)[stim6hMono],paste0("./footprint/1h_vs_6h/",stim,"/",stim,"_6h_Mono.tsv"),sep="\n",row.names=FALSE,col.names = FALSE,quote=FALSE)

# All three conditions
monoCtrlStim <- ATAC.meta$Condition %in% c("Control_1h",paste(stim,c("1h","6h"),sep="_")) & ATAC.meta$pairedLabel2 %in% "Mono"


# Let us define super-responders in the LPS 1hr group, and super slow responders in the LPS 6hr group for mono
# Use NN stim time window cut-off to define this

# NOTE: ORIGINALLY WE USED LOWER CUT 0.3, UPPER CUT 0.7 (scaled) for the DE fig
# But we also tried 0.4 and
# Changed it below to include more cells (useful for footprinting)
lower.cut <- 0.3 # 2 if unscaled
upper.cut <- 0.7 # 4 if unscaled 


SR <- ATAC.meta$Condition %in% paste(stim,"1h",sep="_") & (pseudo.d[,paste0(stim,".scaled")] >= lower.cut & pseudo.d[,paste0(stim,".scaled")] <= upper.cut) & ATAC.meta$pairedLabel2 %in% cellType
table(SR)

SSR <- ATAC.meta$Condition %in% paste(stim,"6h",sep="_") & (pseudo.d[,paste0(stim,".scaled")] >= lower.cut & pseudo.d[,paste0(stim,".scaled")] <= upper.cut) & ATAC.meta$pairedLabel2 %in% cellType
table(SSR)

# NOTE: USED THE ABOVE CELLS TO WRITE BARCODE FILES, AND GET TOTAL PROMOTER COUNTS FOR BIGWIG NORMALIZATION
write.table(rownames(ATAC.meta)[SR],paste0("./footprint/leader_vs_lagger/",stim,"/",stim,"_leaders_",lower.cut,"_",upper.cut,".tsv"),sep="\n",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(rownames(ATAC.meta)[SSR],paste0("./footprint/leader_vs_lagger/",stim,"/",stim,"_laggers_",lower.cut,"_",upper.cut,".tsv"),sep="\n",row.names = FALSE,col.names = FALSE,quote = FALSE)


# Do differential motifs (latest re-run)
dev_motif <- readRDS("./data_freeze_10_16_2020/chromVAR/stim_motif_dev.rds")
dim(dev_motif) # All scATAC cells

dev_motif <- dev_motif[,rownames(ATAC.meta)] # Subset cells, map order

# Don't use bagged here, since we want to know which the real motif activator/repressor might be
Z <- deviationScores(dev_motif)
rownames(Z) <- extractTFNames(rownames(Z))
dim(Z)

# Order motif matrix by nn stim time
# Do this only for Control 1h, stim 1h, and stim 6h monocytes

# Highlight chosen cells
ggplot(ATAC.meta,aes(UMAP1,UMAP2)) + geom_point_rast(color="gray80") + geom_point(data = ATAC.meta[monoCtrlStim,],color="firebrick",size=0.5) + theme_classic()
# Highlight SR and SSR
# Leaders
ggplot(ATAC.meta,aes(UMAP1,UMAP2)) + geom_point_rast(color="gray80") + geom_point(data = ATAC.meta[SR,],color="darkorange",size=0.5) + theme_classic()
# Laggers
ggplot(ATAC.meta,aes(UMAP1,UMAP2)) + geom_point_rast(color="gray80") + geom_point(data = ATAC.meta[SSR,],color="cadetblue",size=0.5) + theme_classic()


# Differential wrt leaders and laggers
#monoBagLeaders <- extractTFNames(names(bagDeviations(dev_motif[,monoCtrlStim],cor = 0.8,organism = "human")))

# Do DE between these two groups
# Pseudo time (should not be different)
t.test(pseudo.d[SR,paste0(stim,".scaled")],pseudo.d[SSR,paste0(stim,".scaled")])
wilcox.test(pseudo.d[SR,paste0(stim,".scaled")],pseudo.d[SSR,paste0(stim,".scaled")])

# Leader vs lagger
DE_motif <- wilcoxDE(Z[,SR],Z[,SSR])

# Also do stim 6h vs stim 1h Mono
DE_motif6h1h <- wilcoxDE(Z[,stim6hMono],Z[,stim1hMono])
DE_motifLPS1hCtrl1h <- wilcoxDE(Z[,stim1hMono],Z[,ctrl1hMono])
all.equal(rownames(DE_motif),rownames(DE_motif6h1h))
all.equal(rownames(DE_motif),rownames(DE_motifCtl1hLPS1h))

motif.d <- data.frame("Motif"=rownames(DE_motif),
                      logFDR_LL=DE_motif$signlogFDR,
                      logFDR_6h1h=DE_motif6h1h$signlogFDR,
                      logFDR_LPS1hCtrl1h=DE_motifLPS1hCtrl1h$signlogFDR,
                      stringsAsFactors = FALSE)

# Highlight as significant if Leader-lagger or 1h 6h is sig at 0.05
isSigMotif <- abs(motif.d$logFDR_LL) >= -log10(0.05) | abs(motif.d$logFDR_6h1h) >= -log10(0.05)
motif.d$sig <- ifelse(isSigMotif,"Yes","No")


TFhighlights <- c("JUNB","JUN","FOSL1","FOS","SMARCC1","NFKB1","NFKB2","CEPBA","SMAD1","SMAD5","NRF1","STAT5A","IRF6","NFYC","STAT3","STAT1","RELA","IRF2","IRF9","ELK1","ELK4","ETV1","GABPA","ETV4","SNAI1","SNAI2","CTCF","IRF3","IRF5","CEBPZ","TCF4","NFYB")
motif.d$Label <- motif.d$Motif
motif.d$Label[!motif.d$Label %in% TFhighlights] <- ""
#motif.d$Label[abs(motif.d$logFDR_6h1h) < 3 & abs(motif.d$logFDR_LL) < -log10(0.05)] <- ""
#motif.d$Label[motif.d$sig=="No"] <- ""

library(ggrepel)
motif.d$dens <- get_density(motif.d$logFDR_6h1h,motif.d$logFDR_LL,100)
gDiffMotif <- ggplot(motif.d,aes(x=logFDR_6h1h,y=logFDR_LL,label=Label,color=dens)) + geom_point(size=0.5,shape=16) + 
  labs(x="log FDR 1h vs 6h",y="log FDR Leader vs Lagger",color="Density")+ 
  scale_color_gradientn(colours = jdb_palette("solar_extra"))+
  geom_text_repel(fontface="italic",size=4,segment.color = "gray",segment.size = 0.5,color="black")+
  #geom_text(fontface="italic",size=3,color="black",nudge_x = 1.1)+
  geom_hline(yintercept = c(-log10(0.05),log10(0.05)),linetype="dashed",color="gray60") + 
  geom_vline(xintercept = c(-log10(0.05),log10(0.05)),linetype="dashed",color="gray60") + 
  theme_classic() + theme(axis.text = element_text(color="black"))


gDiffMotif2 <- ggplot(motif.d,aes(x=logFDR_6h1h,y=logFDR_LL,label=Label,color=sig)) + geom_point(size=0.5,shape=16) + 
  labs(x="log FDR 6h vs 1h",y="log FDR Leader vs Lagger",color="Significant")+ 
  scale_color_manual(values=c("gray60","firebrick"))+
  geom_text_repel(fontface="italic",size=4,segment.color = "gray",segment.size = 0.5,color="black")+
  #geom_text(fontface="italic",size=3,color="black",nudge_x = 1.1)+
  geom_hline(yintercept = c(-log10(0.05),log10(0.05)),linetype="dashed",color="gray60") + 
  geom_vline(xintercept = c(-log10(0.05),log10(0.05)),linetype="dashed",color="gray60") + 
  theme_classic() + theme(axis.text = element_text(color="black"))

pdf("./data_freeze_10_16_2020/figures/scATAC/scATAC_leadervslagger_6hvs1h_Mono_diffMotif_wilcox_unbagged.pdf",height = 5,width = 6,useDingbats = FALSE)
gDiffMotif
gDiffMotif2
dev.off()


# Also see whether some of the top leader vs lagger DE motifs are enriched among Golgi up vs golgi down cells
# Idea being, we want to see if jun-fos activity is variable with golgi response
golgiMonocells <- ATAC.meta$Condition %in% "LPSGolgiPlug_6h" & ATAC.meta$pairedLabel2 %in% cellType
golgiVar <- apply(Z[monoBagLeaders,golgiMonocells],1,sd)
var.d <- data.frame("TF"=names(golgiVar),var=golgiVar) %>% 
  arrange(desc(var)) %>% 
  mutate(TF=factor(TF,levels=as.character(TF))) %>%
  mutate(Label=ifelse(var >= 1.5,as.character(TF),""))

ggplot(var.d,aes(x=TF,y=var,label=Label)) + geom_point(size=0.5) + theme_classic() + 
  theme(axis.text.x = element_blank(),axis.text = element_text(color="black")) + 
  labs(x="Ranked TFs",y="Variability") + geom_text_repel(fontface="italic")

golgi.up <- 0.75
golgi.down <- 0.25

golgi_high_cells <- golgiMonocells & pseudo.d$LPS.scaled >= golgi.up
golgi_low_cells <- golgiMonocells & pseudo.d$LPS.scaled <= golgi.down

write.table(rownames(ATAC.meta)[golgi_high_cells],"./data_freeze_10_16_2020/annot/scATAC_stim_LPS_golgi_high.tsv",sep="\n",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(rownames(ATAC.meta)[golgi_low_cells],"./data_freeze_10_16_2020/annot/scATAC_stim_LPS_golgi_low.tsv",sep="\n",row.names = FALSE,col.names = FALSE,quote = FALSE)

# DE motif golgi high vs golgi low
DE_golgi <- wilcoxDE(Z[,golgi_high_cells],
                     Z[,golgi_low_cells])

gDiffGolgi1 <- DE_golgi %>% as.data.frame() %>% tibble::rownames_to_column(var = "TF") %>% 
  mutate(DE_6h1h=DE_motif6h1h$signlogFDR,DE_LL=DE_motif$signlogFDR) %>%
  mutate(Label=ifelse(abs(DE_LL) >= -log10(0.05) & logFDR >= -log10(0.05),TF,"")) %>%
  mutate(isSig=ifelse(abs(DE_LL) >= -log10(0.05) | logFDR >= -log10(0.05),"Yes","No")) %>%
  ggplot(aes(x=DE_LL,y=signlogFDR,label=Label,color=isSig)) + 
  geom_hline(yintercept = c(-log10(0.05),log10(0.05)),linetype="dashed",color="gray60") + 
  geom_vline(xintercept = c(-log10(0.05),log10(0.05)),linetype="dashed",color="gray60") + 
  geom_point(size=0.5) + theme_classic() + 
  scale_color_manual(values=c("gray60","firebrick"))+
  geom_text_repel(fontface="italic",size=4,segment.color = "gray",segment.size = 0.5,color="black")+
  labs(y="log FDR Golgi high v Golgi low",x="log FDR Leader vs Lagger") + 
  theme(axis.text = element_text(color="black"))
gDiffGolgi1

# Also include DE Golgi with 6h vs 1h
gDiffGolgi2 <- DE_golgi %>% as.data.frame() %>% tibble::rownames_to_column(var = "TF") %>% 
  mutate(DE_6h1h=DE_motif6h1h$signlogFDR,DE_LL=DE_motif$signlogFDR) %>%
  mutate(Label=ifelse((abs(DE_6h1h) >= 3 & logFDR >= -log10(0.05) | TF %in% TFhighlights),TF,"")) %>%
  mutate(isSig=ifelse(abs(DE_6h1h) >= -log10(0.05) | logFDR >= -log10(0.05),"Yes","No")) %>%
  ggplot(aes(x=DE_6h1h,y=signlogFDR,label=Label,color=isSig)) + 
  geom_hline(yintercept = c(-log10(0.05),log10(0.05)),linetype="dashed",color="gray60") + 
  geom_vline(xintercept = c(-log10(0.05),log10(0.05)),linetype="dashed",color="gray60") + 
  geom_point(size=0.5) + theme_classic() + 
  scale_color_manual(values=c("gray60","firebrick"))+
  geom_text_repel(fontface="italic",size=4,segment.color = "gray",segment.size = 0.5,color="black")+
  labs(y="log FDR Golgi high v Golgi low",x="log FDR 6h vs 1h") + 
  theme(axis.text = element_text(color="black"))
gDiffGolgi2


pdf("./data_freeze_10_16_2020/figures/scATAC/scATAC_leadervslagger_6h1h_GolgiHighGolgiLow_Mono_diffMotif_wilcox_unbagged.pdf",height = 5,width = 6,useDingbats = FALSE)
gDiffGolgi1
gDiffGolgi2
dev.off()

