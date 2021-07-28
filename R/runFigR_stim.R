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
dorcMat.smoothed <- smoothGeneScoresNN(lsi.knn,dorcMat,nCores = 6)
stopifnot(all.equal(rownames(dorcMat.smoothed),dorcGenes))

rm(DORC.SE)
gc()

# Just so that the smoothing function will work, since it checks for matching attributes
rownames(lsi.knn) <- colnames(RNA.SE)

# Run only on TFs to save time
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

plotDrivers(figR.d = stim_FigR,marker = "MX1")
plotDrivers(figR.d = stim_FigR,marker = "REL")
plotDrivers(figR.d = stim_FigR,marker = "TRAF1")

# Load DORC-implicated SNP P mat
SNPPmat <- readRDS("./data_freeze_10_16_2020/GWAS/GWAS_DORC_SNPoV_unfiltered_P_mat.rds")
SNPDORCs <- rownames(SNPPmat) 
length(SNPDORCs)

SNPPmat[SNPPmat > 30] <- 30 # Cap

# Heatmap 
figR_heat <- plotfigRHeatmap(figR.d = stim_FigR,
                             score.cut = 1.5,
                             DORCs = SNPDORCs,
                             column_names_gp=gpar(fontsize=5),
                             show_row_dend = FALSE,
                             row_names_side = "left")


myDORCs <- stim_FigR %>% filter(abs(Score)>=1.5 & DORC %in% SNPDORCs) %>% pull(DORC) %>% unique()
myDORCs <- sort(myDORCs) # Has to be the same input as would be for main dorc heatmap (alphabetically sorted since we use reshape)

SNP_heat <- Heatmap(as.matrix(SNPPmat[sort(myDORCs),]),
                    width = 20,
                    border = TRUE,
                    use_raster = FALSE,
                    name = "GWAS log P",
                    show_row_dend = FALSE,
                    cluster_rows = FALSE,
                    show_row_names = TRUE,
                    col = colorRamp2(colors = rev(jdb_palette("china_weirdo",type = "discrete",n = 8))[c(1,5,8)],breaks = c(1,10,20)), # Assuming for log GWAS P
                    row_names_gp = gpar(fontsize = 7,fontface="italic"),
                    column_names_gp = gpar(fontsize = 7),
                    column_names_rot = -30,
                    column_dend_height = unit(0.2,"cm"),
                    clustering_distance_columns = "spearman",
                    heatmap_legend_param = list(direction = "vertical"))
SNP_heat

ht_list <- figR_heat + SNP_heat
draw(ht_list)

# Load DORC list (based on peak Ov) per disease
dorc.list <- readRDS("./data/DORCs/DORCs_Ov_per_disease.rds")


# Just SLE-implicated DORCs
figR_heat_SLE <- plotfigRHeatmap(figR.d = stim_FigR,
                             score.cut = 1.5,
                             DORCs = dorc.list$SystemicLupusErymatosus[!grepl("HLA",dorc.list$SystemicLupusErymatosus)],
                             column_names_gp=gpar(fontsize=5),
                             show_row_dend = FALSE,
                             row_names_side = "left")

# D3Network
plotfigRNetwork(stim_FigR,score.cut = 1.5,DORCs = dorc.list$SystemicLupusErymatosus[!grepl("HLA",dorc.list$SystemicLupusErymatosus)])


# Using ggnet2

lildat <- stim_FigR %>% filter(DORC %in% dorc.list$SystemicLupusErymatosus[!grepl("HLA",dorc.list$SystemicLupusErymatosus)]) %>% filter(abs(Score) >= 1.5)
lildat$Motif <- paste0(lildat$Motif, ".")
lildat$DORC <- paste0(lildat$DORC)
lildat$weight <- scales::rescale(abs(lildat$Score),to=c(0.5,3))

net <- network(lildat[,c("Motif","DORC")],
               directed = TRUE,
               matrix.type="edgelist")

set.edge.value(x = net,attrname = "weight",value = lildat$weight)


# Access node name, and assign as either DORC/Motif
net %v% "class" = ifelse(sapply(net$val,"[[",2) %in% lildat$DORC, "DORC", "Motif")

set.vertex.attribute(net,"vertex.names",gsub(pattern = ".",replacement = "",x = get.vertex.attribute(net,"vertex.names"),fixed = TRUE))

# Assign edge color based on positive or negative correlation
set.edge.attribute(net, "color", ifelse(lildat$Corr > 0, "firebrick", "steelblue"))

# set colors for each mode
col = c("Motif" = "gray60", "DORC" = "sandybrown")

ggnet2(net = net,color.legend = "C",
       size.legend = "Degree",
       label = TRUE,
       size="class",
       #size="outdegree",
       edge.size = "weight",
       edge.color = "color",
       #node.label=NA,
       edge.alpha = 0.75,
       alpha = 0.7,
       color="class",
       size.palette = c("DORC" = 3, "Motif" = 1),
       label.size = 3,
       #palette = "Set2",
       palette=col,
       #arrow.size = 4,
       #arrow.gap = 0.025,
       legend.position = "none") + 
  coord_equal()