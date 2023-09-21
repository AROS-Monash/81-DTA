# WGCNA on Dip selected genes
#' @author Shani Amaraisnghe
#' @date 29/07/2022
#' @aliases Rosello-Diez lab, ARMI, Monash
#' @concept analysing Trev-seq bulk RNA-Seq data for Chee
#' @description for the july-2021 the data was demultiplexed using umi-tools, and bbmap 
#' had to spend two rounds of demultiplexing to get the data seperated to samples and the UMIs.
#' For sept-2021 data the original fastq files available in the raw data directory was used as 
#' they were already demultiplexed
#' for nov-2021 data, one round of umi-tools was carried out to demultiplex by UMI  
#' next, the demux data was aligned read counting was done using featureCounts()
#' next DGEobject was created using edgeR with counts filtered and normalised
#' DE analysis was performed
#' Now I want to look at the reglatory networks of this data and whether I can relate the DE genes to these networks
#' using WGCNA in M3

# ------------- Load libraries -----------
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
library(here)
library(purrr)
library(stringr)
library(edgeR)
library(limma)
library(Rsamtools)
library(ShortRead)
library(Glimma)
library(RColorBrewer)
library(Mus.musculus) 
library(org.Mm.eg.db)
library(WGCNA)
library(diptest)
library(viridis)
library(gplots)

# ------------------ E15.5 ---------------------------------------------
# ----------------------------------------------------------------------

# ------------- Set working directory ----
# location of processing
setwd("/fs03/ag36/Shani/Trev-Seq/featurecounts_edgeR_Mar_2023_dedup/E15")
here::here()
dir.create("./WGCNA")
dir.create("./WGCNA/E15-5") # this was needed because I have set up the helper functions to have this format
source("/fs04/ag36/Shani/Trev-Seq/scripts/helper_scripts.R")
source("/fs04/ag36/Shani/Trev-Seq/scripts/WGCNA_helper_functions.R")

metadata <- read.table("/fs04/ag36/Shani/Trev-Seq/fastq_exp_dedup.csv", sep = "\t", header = TRUE)

# getting some more physiological data in the metadata 
metadata$GENO <- factor(paste0(substr(metadata$GENOTYPE, 1,1), 
                        "_",
                        substr(metadata$LIMB, 1,1)))


# --------------- Load data -----------------
dge.15    <- readRDS("./dge15.final.rds")
v.15      <- readRDS("./voom15.final.rds")
efit.15   <- readRDS("./efit15.final.rds")
design.15 <- readRDS("./design15.final.rds")
tt.15     <- read.csv("./toptable_15_p1.csv")

# ----------- run diptest ----------------
dip <- dip.test(dge.15$counts)

# > (dip)

#         Hartigans' dip test for unimodality / multimodality

# data:  dge.15$counts
# D = 0.0050998, p-value < 2.2e-16
# alternative hypothesis: non-unimodal, i.e., at least bimodal

#------------- subset tt.15 and setup inputs --------------------

tt.15     <- tt.15[tt.15$adj.P.Val <= 2.2e-16,]
counts.15 <- cpm(dge.15$counts[row.names(tt.15),])
traits.15 <- metadata$expgroup[match(colnames(counts.15), metadata$ID.NEW)]
stage.15  <- "E15-5"

# ------------------ nthreads ------------------
Nthreads <- 12
enableWGCNAThreads(nthread = Nthreads)

setwd("/fs03/ag36/Shani/Trev-Seq/featurecounts_edgeR_2023_dedup/E15/WGCNA/") # this was needed because I have set up the helper functions to have this format

y.15           <- load_and_define(counts.15, stage.15, traits.15, metadata)

euc.dist.plot(y.15, stage.15, traits.15)
powers.15      <- sft_threshold(y.15, stage.15)
# calculate power.15 by non cpm matrix
# soft threshold selected: 9
sft.15         <- 9
k.15           <- k_calc(y.15, sft.15, stage.15) 
#   scaleFreeRsquared slope
# 1              0.91 -1.38
net.15         <- wgcna.automatic(y.15, sft.15, stage.15, 0.95)
tom.15         <- wgcna.tom(y.15, sft.15, stage.15, 12)
disstom.15     <- wgcna.disstom(tom.15, stage.15, sft.15)
merged.15      <- construct.modules(disstom.15, y.15, stage.15, 30, traits.15, tom.15, sft.15, 0.5, 0.6)

ME_heatmaps(merged.15, stage.15)

# ------------------ E17.5 ---------------------------------------------
# ----------------------------------------------------------------------

# ------------- Set working directory ----
# location of processing
setwd("/fs03/ag36/Shani/Trev-Seq/featurecounts_edgeR_2023_dedup/E17")
here::here()
dir.create("./WGCNA")
dir.create("./WGCNA/E17-5") # this was needed because I have set up the helper functions to have this format
source("/fs04/ag36/Shani/Trev-Seq/scripts/helper_scripts.R")
source("/fs04/ag36/Shani/Trev-Seq/scripts/WGCNA_helper_functions.R")

metadata <- read.table("/fs04/ag36/Shani/Trev-Seq/fastq_exp_dedup.csv", sep = "\t", header = TRUE)

# getting some more physiological data in the metadata 
metadata$GENO <- factor(paste0(substr(metadata$GENOTYPE, 1,1), 
                        "_",
                        substr(metadata$LIMB, 1,1)))


# --------------- Load data -----------------
dge.17    <- readRDS("./dge17.final.rds")
v.17      <- readRDS("./voom17.final.rds")
efit.17   <- readRDS("./efit17.final.rds")
design.17 <- readRDS("./design17.final.rds")
tt.17     <- read.csv("./toptable_17_CC_new_2.csv")

# ----------- run diptest ----------------
dip <- dip.test(dge.17$counts)

# > dip

#         Hartigans' dip test for unimodality / multimodality

# data:  dge.17$counts
# D = 0.0080067, p-value < 2.2e-16
# alternative hypothesis: non-unimodal, i.e., at least bimodal

#------------- subset tt.17 and setup inputs --------------------

tt.17     <- tt.17[tt.17$adj.P.Val <= 2.2e-16,]
counts.17 <- (dge.17$counts[row.names(tt.17),])
traits.17 <- metadata$expgroup[match(colnames(counts.17), metadata$ID.NEW)]
stage.17  <- "E17-5"

# ------------------ nthreads ------------------
Nthreads <- 12
enableWGCNAThreads(nthread = Nthreads)

setwd("/fs03/ag36/Shani/Trev-Seq/featurecounts_edgeR_2023_dedup/E17/WGCNA/") # this was needed because I have set up the helper functions to have this format

y.17           <- load_and_define(counts.17, stage.17, traits.17, metadata)

euc.dist.plot(y.17, stage.17, traits.17)
powers.17      <- sft_threshold(y.17, stage.17)
# calculate power.17 by non cpm matrix

# soft threshold selected: 9
sft.17         <- 9
k.17           <- k_calc(y.17, sft.17, stage.17) 
#   scaleFreeRsquared slope
# 1              0.91 -1.38

counts.17      <- cpm(dge.17$counts[row.names(tt.17),])
net.17         <- wgcna.automatic(y.17, sft.17, stage.17, 0.95)
tom.17         <- wgcna.tom(y.17, sft.17, stage.17, 12)
disstom.17     <- wgcna.disstom(tom.17, stage.17, sft.17)
merged.17      <- construct.modules(disstom.17, y.17, stage.17, 30, traits.17, tom.17, sft.17, 0.5, 0.6)

ME_heatmaps(merged.17, stage.17)


