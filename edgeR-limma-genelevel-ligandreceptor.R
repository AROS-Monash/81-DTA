#' ----- Gene set testing on Trev-Seq dataset (gene level) --------
#' @author Shani Amaraisnghe
#' @date 09/08/2022
#' @aliases Rosello-Diez lab, ARMI, Monash
#' @concept analysing Trev-seq bulk RNA-Seq data for Chee
#' @description for the july-2021 the data was demultiplexed using umi-tools, and bbmap 
#' had to spend two rounds of demultiplexing to get the data seperated to samples and the UMIs.
#' For sept-2021 data the original fastq files available in the raw data directory was used as 
#' they were already demultiplexed
#' for nov-2021 data, one round of umi-tools was carried out to demultiplex by UMI  
#' next, the demux data was aligned and read counting was done using featureCounts()
#' then edgeR-limma pipeline was run on gene-level
#' now checking for ligand-receptor interaction related genes in the DE lists
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
library(EDASeq)
library(Glimma)
library(RColorBrewer)
library(Mus.musculus) 
library(org.Mm.eg.db)

# -------------I/O --------------
# working in M3
#module load R/4.2-mkl

# ------------- Set working directory ----

# location of processing
setwd("/fs03/ag36/Shani/Trev-Seq/featurecounts_edgeR")
here()
dir.create("./LRpairs")

# --------------- Load data -----------------
dge.15 <- readRDS("dge.15.final.rds")
dge.17 <- readRDS("dge.17.final.rds")
dge.0 <- readRDS("dge.0.final.rds")
dge.3 <- readRDS("dge.3.final.rds")

metadata <- read.table("/fs02/ag36/Shani/Trev-Seq/fastq_exp.csv", sep = ",", header = TRUE)

metadata$expgroup <- paste(substr(metadata$GENOTYPE, 1,1),
                            metadata$SIDE,
                            substr(metadata$STAGE, 2,3),
                            sep = "_")

metadata$ID.NEW  <- factor(paste(metadata$MONTH, metadata$SAMPLE, sep ="_"))

metadata$ID.2    <- factor(paste(metadata$MONTH, metadata$SAMPLE, metadata$REPLICATE, sep ="_"))

# -------------- Load ligand-receptor pair information ------------

# downloaded from http://tcm.zju.edu.cn/celltalkdb/download.php

lrpairs <- read.table((file.path("LRpairs", "mouse_lr_pair.txt")), header = TRUE, sep ="\t")

# load the toptables
toptable.15 <- read.table("./toptable_15_p1.csv", sep=",", header= TRUE)
toptable.17 <- read.table("./toptable_17_p1.csv", sep=",", header= TRUE)
toptable.0  <- read.table("./toptable_0_p1.csv", sep=",", header= TRUE)
toptable.3  <- read.table("./toptable_3_p1.csv", sep=",", header= TRUE)

# find the common ligands - 15.5
common.ligand.15     <- subset(lrpairs, lrpairs$ligand_ensembl_gene_id %in% toptable.15$efit.15.coef.match.row.names.out.15.)
common.receptor.15   <- subset(lrpairs, lrpairs$receptor_ensembl_gene_id %in% toptable.15$efit.15.coef.match.row.names.out.15.)

write.table(common.ligand.15, "./LRpairs/common.ligand.15", sep=",", quote = FALSE, row.names = TRUE)
write.table(common.receptor.15, "./LRpairs/common.receptor.15", sep=",", quote = FALSE, row.names = TRUE)

# find the common ligands - 17.5
common.ligand.17     <- subset(lrpairs, lrpairs$ligand_ensembl_gene_id %in% toptable.17$efit.17.coef.match.row.names.out.17.)
common.receptor.17   <- subset(lrpairs, lrpairs$receptor_ensembl_gene_id %in% toptable.17$efit.17.coef.match.row.names.out.17.)

write.table(common.ligand.17, "./LRpairs/common.ligand.17", sep=",", quote = FALSE, row.names = TRUE)
write.table(common.receptor.17, "./LRpairs/common.receptor.17", sep=",", quote = FALSE, row.names = TRUE)

# find the common ligands - P0
common.ligand.0     <- subset(lrpairs, lrpairs$ligand_ensembl_gene_id %in% toptable.0$efit.0.coef.match.row.names.out.0.)
common.receptor.0   <- subset(lrpairs, lrpairs$receptor_ensembl_gene_id %in% toptable.0$efit.0.coef.match.row.names.out.0.)

write.table(common.ligand.0, "./LRpairs/common.ligand.0", sep=",", quote = FALSE, row.names = TRUE)
write.table(common.receptor.0, "./LRpairs/common.receptor.0", sep=",", quote = FALSE, row.names = TRUE)

# find the common ligands - P3
common.ligand.3     <- subset(lrpairs, lrpairs$ligand_ensembl_gene_id %in% toptable.3$efit.3.coef.match.row.names.out.3.)
common.receptor.3   <- subset(lrpairs, lrpairs$receptor_ensembl_gene_id %in% toptable.3$efit.3.coef.match.row.names.out.3.)

write.table(common.ligand.3, "./LRpairs/common.ligand.3", sep=",", quote = FALSE, row.names = TRUE)
write.table(common.receptor.3, "./LRpairs/common.receptor.3", sep=",", quote = FALSE, row.names = TRUE)

# ------------- for the CC results ------------------------------------

# load the toptables
toptable.15 <- read.table("./toptable_15_Nov_p1.CC.csv", sep=",", header= TRUE)
toptable.17 <- read.table("./toptable_17_p1.CC.csv", sep=",", header= TRUE)
toptable.0  <- read.table("./toptable_0_p1.CC.csv", sep=",", header= TRUE)
toptable.3  <- read.table("./toptable_3_p1.CC.csv", sep=",", header= TRUE)

# find the common ligands - 15.5
common.ligand.15     <- subset(lrpairs, lrpairs$ligand_ensembl_gene_id %in% toptable.15$efit.15.coef.match.row.names.out.15.)
common.ligand.15     <- merge(common.ligand.15, toptable.15,
                                by.x = "ligand_ensembl_gene_id",
                                by.y = "efit.15.coef.match.row.names.out.15.", 
                                all.x = TRUE)

common.receptor.15   <- subset(lrpairs, lrpairs$receptor_ensembl_gene_id %in% toptable.15$efit.15.coef.match.row.names.out.15.)
common.receptor.15   <- merge(common.receptor.15, toptable.15,
                                by.x = "receptor_ensembl_gene_id",
                                by.y = "efit.15.coef.match.row.names.out.15.", 
                                all.x = TRUE)

write.table(common.ligand.15, "./LRpairs/common.ligand.15.CC.tab", sep="\t", quote = FALSE, row.names = FALSE)
write.table(common.receptor.15, "./LRpairs/common.receptor.15.CC.tab", sep="\t", quote = FALSE, row.names = FALSE)

# find the common ligands - 17.5
common.ligand.17     <- subset(lrpairs, lrpairs$ligand_ensembl_gene_id %in% toptable.17$efit.17.coef.match.row.names.out.17.)
common.ligand.17     <- merge(common.ligand.17, toptable.17,
                                by.x = "ligand_ensembl_gene_id",
                                by.y = "efit.17.coef.match.row.names.out.17.", 
                                all.x = TRUE)

common.receptor.17   <- subset(lrpairs, lrpairs$receptor_ensembl_gene_id %in% toptable.17$efit.17.coef.match.row.names.out.17.)
common.receptor.17   <- merge(common.receptor.17, toptable.17,
                                by.x = "receptor_ensembl_gene_id",
                                by.y = "efit.17.coef.match.row.names.out.17.", 
                                all.x = TRUE)

write.table(common.ligand.17, "./LRpairs/common.ligand.17.CC.tab", sep="\t", quote = FALSE, row.names = FALSE)
write.table(common.receptor.17, "./LRpairs/common.receptor.17.CC.tab", sep="\t", quote = FALSE, row.names = FALSE)

# find the common ligands - P0
common.ligand.0     <- subset(lrpairs, lrpairs$ligand_ensembl_gene_id %in% toptable.0$efit.0.coef.match.row.names.out.0.)
common.ligand.0     <- merge(common.ligand.0, toptable.0,
                                by.x = "ligand_ensembl_gene_id",
                                by.y = "efit.0.coef.match.row.names.out.0.", 
                                all.x = TRUE)

common.receptor.0   <- subset(lrpairs, lrpairs$receptor_ensembl_gene_id %in% toptable.0$efit.0.coef.match.row.names.out.0.)
common.receptor.0   <- merge(common.receptor.0, toptable.0,
                                by.x = "receptor_ensembl_gene_id",
                                by.y = "efit.0.coef.match.row.names.out.0.", 
                                all.x = TRUE)

write.table(common.ligand.0, "./LRpairs/common.ligand.0.CC.tab", sep="\t", quote = FALSE, row.names = FALSE)
write.table(common.receptor.0, "./LRpairs/common.receptor.0.CC.tab", sep="\t", quote = FALSE, row.names = FALSE)

# find the common ligands - P3
common.ligand.3     <- subset(lrpairs, lrpairs$ligand_ensembl_gene_id %in% toptable.3$efit.3.coef.match.row.names.out.3.)
common.ligand.3     <- merge(common.ligand.3, toptable.3,
                                by.x = "ligand_ensembl_gene_id",
                                by.y = "efit.3.coef.match.row.names.out.3.", 
                                all.x = TRUE)

common.receptor.3   <- subset(lrpairs, lrpairs$receptor_ensembl_gene_id %in% toptable.3$efit.3.coef.match.row.names.out.3.)
common.receptor.3   <- merge(common.receptor.3, toptable.3,
                                by.x = "receptor_ensembl_gene_id",
                                by.y = "efit.3.coef.match.row.names.out.3.", 
                                all.x = TRUE)

write.table(common.ligand.3, "./LRpairs/common.ligand.3.CC.tab", sep="\t", quote = FALSE, row.names = FALSE)
write.table(common.receptor.3, "./LRpairs/common.receptor.3.CC.tab", sep="\t", quote = FALSE, row.names = FALSE)
