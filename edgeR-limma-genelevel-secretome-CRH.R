
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
#' now checking for secretome related genes in the DE lists
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
library(Glimma)
library(RColorBrewer)
library(Mus.musculus) 
library(org.Mm.eg.db)

# -------------I/O --------------
# working in M3
#module load R/4.2-mkl

# ------------- Set working directory ----

# location of processing
setwd("/fs03/ag36/Shani/Trev-Seq/CRH")
here::here()
dir.create("./mouse_secretome")

# --------------- Load data -----------------

crh.top.24 <- read.table(file.path("/fs03/ag36/Shani/Trev-Seq", "CRH", "toptable_Unx.24vsSham.24_crh.csv"), sep=",", header= TRUE)
crh.top.48 <- read.table(file.path("/fs03/ag36/Shani/Trev-Seq", "CRH", "toptable_Unx.48vsSham.48_crh.csv"), sep=",", header= TRUE)
crh.top.72 <- read.table(file.path("/fs03/ag36/Shani/Trev-Seq", "CRH", "toptable_Unx.72vsSham.72_crh.csv"), sep=",", header= TRUE)

# -------------- Load mouse secretome information ------------

mouse_secretome<- read.table((file.path("/fs03/ag36/Shani/Trev-Seq", "featurecounts_edgeR", "mouse_secretome", "secretome_results.csv")), header = TRUE, sep =",")

# ------------- Findng the common secretome genes/protiens ---------------------

# find - 24 HR
mouse.secretome.24     <- subset(mouse_secretome, mouse_secretome$Gene %in% crh.top.24$genes)
write.table(mouse.secretome.24, "./mouse_secretome/mouse.secretome.24", sep=",", quote = FALSE, row.names = TRUE)

# find - 48 HR
mouse.secretome.48     <- subset(mouse_secretome, mouse_secretome$Gene %in% crh.top.48$genes)
write.table(mouse.secretome.48, "./mouse_secretome/mouse.secretome.48", sep=",", quote = FALSE, row.names = TRUE)

# find - 72 HR
mouse.secretome.72     <- subset(mouse_secretome, mouse_secretome$Gene %in% crh.top.72$genes)
write.table(mouse.secretome.72, "./mouse_secretome/mouse.secretome.72", sep=",", quote = FALSE, row.names = TRUE)
