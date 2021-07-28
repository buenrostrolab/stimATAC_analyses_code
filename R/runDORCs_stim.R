# Script to run core DORC calling function on stimulation data

### Author: Vinay Kartha
### Contact: <vinay_kartha@g.harvard.edu>
### Affiliation: Buenrostro Lab, Department of Stem Cell and Regenerative Biology, Harvard University

setwd("<data_analysis_folder>")


source("./code/DORC_functions.R")

# ---------------------------------------------------------------------------------- Load stim data
# Load SE object
SE.filt <- readRDS("./data/SE/ATAC_stim_paired.rds")
SE.filt


# Normalized counts from SEs
# Load paired RNA
RNA.SE <- readRDS("./data/SE/RNA_stim_paired.rds")
names(assays(RNA.SE))
rnaMat <- assays(RNA.SE)$counts.norm


# ---------------------------------------------------------------------------------- Call gene-peak associations

# Use 8 cores here
cisCor <- runGenePeakcorr(ATAC.se = SE.filt,
                           RNAmat = rnaMat,
                           genome = "hg19",
                          windowPadSize = 50000,
                           nCores = 8,
                           n_bg = 100,
                           p.cut = NULL)
 
# Filter associations using correlation p-value 
cisCor.filt <- cisCor %>% filter(pvalZ <= 0.05)


# ---------------------------------------------------------------------------------- Determine DORC genes

# Make J plot to rank genes by # significant peaks
dorcGenes <- dorcJplot(dorcTab = cisCor.filt,
                       cutoff = 7,
                       returnGeneList = TRUE)