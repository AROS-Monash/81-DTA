#' ----- Gene set testing on Trev-Seq dataset (gene level) --------
#' @author Shani Amaraisnghe
#' @date 15/08/2022
#' @aliases Rosello-Diez lab, ARMI, Monash
#' @concept analysing CRH bulk RNA-Seq data for comparison with TrevSeq data
#' @description for the july-2021 the data was demultiplexed using umi-tools, and bbmap 
#' had to spend two rounds of demultiplexing to get the data seperated to samples and the UMIs.
#' For sept-2021 data the original fastq files available in the raw data directory was used as 
#' they were already demultiplexed
#' for nov-2021 data, one round of umi-tools was carried out to demultiplex by UMI  
#' next, the demux data was aligned and read counting was done using featureCounts()
#' then edgeR-limma pipeline was run on gene-level
#' now checking for ligand-receptor interaction related genes in the DE lists of CRH dataset
#' 

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
library(dplyr)
library(Glimma)
library(RColorBrewer)
library(Mus.musculus) 
library(org.Mm.eg.db)
library(tidyverse, warn.conflicts = FALSE)

# -------------I/O --------------
# working in M3
#module load R/4.2-mkl

# ------------- Set working directory ----

# location of processing
setwd("/fs03/ag36/Shani/Trev-Seq/CRH")
here::here()
dir.create("./LRpairs")

# --------------- Load data -----------------

crh.top.24 <- read.table(file.path("/fs03/ag36/Shani/Trev-Seq", "CRH", "toptable_Unx.24vsSham.24_crh.csv"), sep=",", header= TRUE)
crh.top.48 <- read.table(file.path("/fs03/ag36/Shani/Trev-Seq", "CRH", "toptable_Unx.48vsSham.48_crh.csv"), sep=",", header= TRUE)
crh.top.72 <- read.table(file.path("/fs03/ag36/Shani/Trev-Seq", "CRH", "toptable_Unx.72vsSham.72_crh.csv"), sep=",", header= TRUE)

# downloaded from http://tcm.zju.edu.cn/celltalkdb/download.php

lrpairs <- read.table((file.path("/fs03/ag36/Shani/Trev-Seq", "featurecounts_edgeR", "LRpairs", "mouse_lr_pair.txt")), header = TRUE, sep ="\t")

# find the common ligands - 24 HR
common.ligand.24     <- subset(lrpairs, lrpairs$ligand_ensembl_gene_id %in% crh.top.24$genes)
common.receptor.24   <- subset(lrpairs, lrpairs$receptor_ensembl_gene_id %in% crh.top.24$genes)

write.table(common.ligand.24, "./LRpairs/common.ligand.24", sep=",", quote = FALSE, row.names = TRUE)
write.table(common.receptor.24, "./LRpairs/common.receptor.24", sep=",", quote = FALSE, row.names = TRUE)

# find the common ligands - 48 HR
common.ligand.48     <- subset(lrpairs, lrpairs$ligand_ensembl_gene_id %in% crh.top.48$genes)
common.receptor.48   <- subset(lrpairs, lrpairs$receptor_ensembl_gene_id %in% crh.top.48$genes)

write.table(common.ligand.48, "./LRpairs/common.ligand.48", sep=",", quote = FALSE, row.names = TRUE)
write.table(common.receptor.48, "./LRpairs/common.receptor.48", sep=",", quote = FALSE, row.names = TRUE)

# find the common ligands - 72 HR
common.ligand.72     <- subset(lrpairs, lrpairs$ligand_ensembl_gene_id %in% crh.top.72$genes)
common.receptor.72   <- subset(lrpairs, lrpairs$receptor_ensembl_gene_id %in% crh.top.72$genes)

write.table(common.ligand.72, "./LRpairs/common.ligand.72", sep=",", quote = FALSE, row.names = TRUE)
write.table(common.receptor.72, "./LRpairs/common.receptor.72", sep=",", quote = FALSE, row.names = TRUE)