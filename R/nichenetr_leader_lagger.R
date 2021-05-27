library(CellChat)
library(nichenetr)
library(Seurat)
library(SummarizedExperiment)
library(ggplot2)
library(dplyr)

# Preprocess data
# Load paired RNA
RNA.SE <- readRDS("/n/holylfs/LABS/buenrostro_lab/Users/vkartha/stim/data_freeze_10_16_2020/SE/RNA_stim_paired.rds")

# Load paired ATAC
ATAC.SE <- readRDS("/n/holylfs/LABS/buenrostro_lab/Users/vkartha/stim/data_freeze_10_16_2020/SE/ATAC_stim_paired.rds")

stopifnot(ncol(ATAC.SE)==ncol(RNA.SE))

ATAC.meta <- as.data.frame(colData(ATAC.SE),stringsAsFactors=FALSE)

nnTime <- read.table("/n/holylfs/LABS/buenrostro_lab/Users/vkartha/stim/data_freeze_10_16_2020/nnTime/LSI_NN_stim_pseudotimev2.tsv",sep="\t",header=TRUE,stringsAsFactors = FALSE)

leader_cells <- rownames(ATAC.meta)[ATAC.meta$Condition %in% "LPS_1h" & ATAC.meta$pairedLabel2 %in% "Mono" & nnTime$LPS.scaled >= 0.3]
lagger_cells <- rownames(ATAC.meta)[ATAC.meta$Condition %in% c("LPS_1h","Control_1h") & ATAC.meta$pairedLabel2 %in% "Mono" & nnTime$LPS.scaled < 0.1]
length(leader_cells)
length(lagger_cells)

# Load leader and lagger cell groups
# This was copied over from the workstation (see /mnt/Analyses/biorad/stim_06_17_2019/footprint/leader_vs_lagger/LPS)

# LPS 1h leader cells
#leader_cells <- read.table("/n/holylfs/LABS/buenrostro_lab/Users/vkartha/stim/footprinting/leaders_vs_laggers/LPS/LPS_leaders_0.3_0.7.tsv",sep="\n",header=FALSE,stringsAsFactors=FALSE)$V1
# LPS 6h lagger cells
#lagger_cells <- read.table("/n/holylfs/LABS/buenrostro_lab/Users/vkartha/stim/footprinting/leaders_vs_laggers/LPS/LPS_laggers_0.3_0.7.tsv",sep="\n",stringsAsFactors=FALSE,header=FALSE)$V1

RNA.mat <- assays(RNA.SE)$counts.norm
dim(RNA.mat)
RNA.mat <- RNA.mat[Matrix::rowSums(RNA.mat)!=0,]

leader.mat <- RNA.mat[,which(rownames(ATAC.meta) %in% leader_cells)]
lagger.mat <- RNA.mat[,which(rownames(ATAC.meta) %in% lagger_cells)]

dim(leader.mat)
dim(lagger.mat)

colnames(leader.mat) <- leader_cells
colnames(lagger.mat) <- lagger_cells

group <- factor(rep(c("leader","lagger"),c(length(leader_cells),length(lagger_cells))))
table(group)

rna.meta <- data.frame(group,row.names = c(leader_cells,lagger_cells))

RNA.seu <- CreateSeuratObject(counts = cbind(leader.mat,lagger.mat),assay = "RNA",meta.data = rna.meta,min.cells = 1)
head(RNA.seu@meta.data)

ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
ligand_target_matrix[1:5,1:5] # target genes in rows, ligands in columns


lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
head(lr_network)


weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))

head(weighted_networks$lr_sig) # interactions and their weights in the ligand-receptor + signaling network
head(weighted_networks$gr) # interactions and their weights in the gene regulatory network

# Define a “sender/niche” cell population and a “receiver/target” cell population present in your expression data and determine which genes are expressed in both populations

# We have to change the default ident to be leader or lagger group
table(RNA.seu$orig.ident)
RNA.seu <- SetIdent(RNA.seu,value = "group")

RNA.seu %>% Idents() %>% unique()

## receiver
receiver = "lagger"
expressed_genes_receiver = get_expressed_genes(receiver, RNA.seu, pct = 0.05)

background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

## sender
sender_celltypes = "leader"

list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, RNA.seu, 0.05) # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()

# Define a gene set of interest: these are the genes in the “receiver/target” cell population that are potentially affected by ligands expressed by interacting cells (e.g. genes differentially expressed upon cell-cell interaction)

# Define gene set of interest
DE_table_receiver = FindMarkers(object = RNA.seu, 
                                ident.1 = "leader", 
                                ident.2 = "lagger", 
                                min.pct = 0.05) %>% tibble::rownames_to_column("gene")

geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.01 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
#geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05) %>% pull(gene)
geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
length(geneset_oi)


RNA.seu <- Seurat::ScaleData(RNA.seu)
Seurat::DoHeatmap(RNA.seu,group.by = "group",features = geneset_oi,assay = "RNA")

# Define a set of potential ligands: these are ligands that are expressed by the “sender/niche” cell population and bind a (putative) receptor expressed by the “receiver/target” population

ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()

expressed_ligands = intersect(ligands,expressed_genes_sender)
expressed_receptors = intersect(receptors,expressed_genes_receiver)

# Perform NicheNet ligand activity analysis
# Rank the potential ligands based on the presence of their target genes in the gene set of interest (compared to the background set of genes)
potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()

ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)

ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))
ligand_activities


# Visualize top 20/25 upstream ligands
best_upstream_ligands = ligand_activities %>% top_n(25, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()

DotPlot(RNA.seu, features = best_upstream_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis()

set.seed(123)
topLigandHeat <- Seurat::DoHeatmap(RNA.seu,group.by = "group",features = best_upstream_ligands,assay = "RNA",
                                   lines.width = 10,raster = FALSE,group.colors = c("gray","firebrick3"),cells = c(sample(lagger_cells,100),leader_cells))
topLigandHeat

# Infer receptors and top-predicted target genes of ligands that are top-ranked in the ligand activity analysis
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 300) %>% bind_rows() %>% tidyr::drop_na()

active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.33)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) + 
  #scale_fill_gradient2(low = "whitesmoke",  high = "darkorange", breaks = c(0,0.0045,0.0090))
  scale_fill_gradient2(low = jdb_palette("brewer_violet")[1],high=jdb_palette("brewer_violet")[9], breaks = c(0,0.0045,0.0090))
p_ligand_target_network

# Top receptors of top-ranked ligands
lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

lr_network_top_df = lr_network_top_df_large %>% tidyr::spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))

vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()

p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")
p_ligand_receptor_network

# Combine ligand activities with DE information
ligand_activities_de = ligand_activities %>% select(test_ligand, pearson) %>% dplyr::rename(ligand = test_ligand) %>% left_join(DE_table_receiver %>% dplyr::rename(ligand = gene))
ligand_activities_de[is.na(ligand_activities_de)] = 0

# make LFC heatmap
lfc_matrix = ligand_activities_de  %>% select(-ligand, -pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities_de$ligand)
rownames(lfc_matrix) = rownames(lfc_matrix) %>% make.names()

order_ligands = order_ligands[order_ligands %in% rownames(lfc_matrix)]
vis_ligand_lfc = lfc_matrix[order_ligands,]

colnames(vis_ligand_lfc) = vis_ligand_lfc %>% colnames() %>% make.names()

p_ligand_lfc = vis_ligand_lfc[,c("p_val","avg_log2FC")] %>% make_threecolor_heatmap_ggplot("Prioritized ligands","LFC in Sender", low_color = "midnightblue",mid_color = "white", mid = median(vis_ligand_lfc), high_color = "red",legend_position = "top", x_axis_position = "top", legend_title = "LFC") + theme(axis.text.y = element_text(face = "italic"))
p_ligand_lfc


# Summarize visualiations

# ligand activity heatmap
ligand_pearson_matrix = ligand_activities %>% select(pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)

rownames(ligand_pearson_matrix) = rownames(ligand_pearson_matrix) %>% make.names()
colnames(ligand_pearson_matrix) = colnames(ligand_pearson_matrix) %>% make.names()

vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")
p_ligand_pearson = vis_ligand_pearson %>% make_heatmap_ggplot("Prioritized ligands","Ligand activity", color = "darkorange",legend_position = "top", x_axis_position = "top", legend_title = "Pearson correlation coefficient\ntarget gene prediction ability)") + theme(legend.text = element_text(size = 9))


# ligand expression Seurat dotplot
order_ligands_adapted = order_ligands
order_ligands_adapted[order_ligands_adapted == "HLA.DRA"] = "HLA-DRA" # cf required use of make.names for heatmap visualization | this is not necessary if these ligands are not in the list of prioritized ligands!
rotated_dotplot = DotPlot(RNA.seu %>% subset(group %in% sender_celltypes), features = order_ligands_adapted, cols = "RdYlBu") + coord_flip() + theme(legend.text = element_text(size = 10), legend.title = element_text(size = 12)) # flip of coordinates necessary because we want to show ligands in the rows when combining all plots

figures_without_legend = cowplot::plot_grid(
  p_ligand_pearson + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()),
  rotated_dotplot + theme(legend.position = "none", axis.ticks = element_blank(), axis.title.x = element_text(size = 12), axis.text.y = element_text(face = "italic", size = 9), axis.text.x = element_text(size = 9,  angle = 90,hjust = 0)) + ylab("Expression in Sender") + xlab("") + scale_y_discrete(position = "right"),
  p_ligand_lfc + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()) + ylab(""),
  p_ligand_target_network + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
  align = "hv",
  nrow = 1,
  rel_widths = c(ncol(vis_ligand_pearson)+6, ncol(vis_ligand_lfc) + 7, ncol(vis_ligand_lfc) + 8, ncol(vis_ligand_target)))

legends = cowplot::plot_grid(
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_pearson)),
  ggpubr::as_ggplot(ggpubr::get_legend(rotated_dotplot)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_lfc)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target_network)),
  nrow = 1,
  align = "h", rel_widths = c(1.5, 1, 1, 1))

combined_plot = cowplot::plot_grid(figures_without_legend, legends, rel_heights = c(10,5), nrow = 2, align = "hv")
combined_plot

# Save plots
pdf("/n/holylfs/LABS/buenrostro_lab/Users/vkartha/stim/data_freeze_10_16_2020/figures/NicheNet/top_ligand_target_network.pdf",width = 10,height = 5)
p_ligand_target_network
dev.off()

pdf("/n/holylfs/LABS/buenrostro_lab/Users/vkartha/stim/data_freeze_10_16_2020/figures/NicheNet/top_ligand_leader_lagger_heat.pdf",width = 8,height = 4)
topLigandHeat
dev.off()
