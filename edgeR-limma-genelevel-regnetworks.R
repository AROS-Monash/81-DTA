#' ----- Gene set testing on Trev-Seq dataset (gene level) --------
#' @author Shani Amaraisnghe
#' @date 08/08/2022
#' @aliases Rosello-Diez lab, ARMI, Monash
#' @concept analysing Trev-seq bulk RNA-Seq data for Chee
#' @description for the july-2021 the data was demultiplexed using umi-tools, and bbmap 
#' had to spend two rounds of demultiplexing to get the data seperated to samples and the UMIs.
#' For sept-2021 data the original fastq files available in the raw data directory was used as 
#' they were already demultiplexed
#' for nov-2021 data, one round of umi-tools was carried out to demultiplex by UMI  
#' next, the demux data was aligned and read counting was done using featureCounts()
#' then edgeR-limma pipeline was run on gene-level
#' running regulatory network analysis using =RTN package for further downstream TFBS analysis
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
library(RTN)

# -------------I/O --------------
# working in M3
#module load R/4.2-mkl

# ------------- Set working directory ----

# location of processing
setwd("/fs03/ag36/Shani/Trev-Seq/featurecounts_edgeR")
here()
dir.create("./Regulons")

# -------------- Load DGE objects ------------
dge.15 <- readRDS("dge.15.Nov.final.rds")
dge.17 <- readRDS("dge.17.final.rds")
dge.0 <- readRDS("dge.0.final.rds")
dge.3 <- readRDS("dge.3.final.rds")

metadata <- readRDS("./final.metadata.rds")

# load the toptables
toptable.15 <- read.table("./toptable_15_Nov_p1.csv", sep=",", header= TRUE)
toptable.17 <- read.table("./toptable_17_p1.csv", sep=",", header= TRUE)
toptable.0  <- read.table("./toptable_0_p1.csv", sep=",", header= TRUE)
toptable.3  <- read.table("./toptable_3_p1.csv", sep=",", header= TRUE)

# ----------- Set up specific data structures -------------

Pappa_TFs <- read.table("/fs02/ag36/Shani/Trev-Seq/PAPPA2_regulators_edited.tab", header = FALSE)
tfs       <- c(Pappa_TFs$V1)

# here TFs should be within the dataset we enter, because in the code they are looking to see whether they are there
# so I can use for example EGR1 and maybe ligands as TFs?
# so I input the protein secuences of the DE genes to the AnimalTFBS database and 
# looked for the TFs there
# dowloaded the TFs list and compared the proteins with the proteins from each DE list 

TF_FULL <- read.table("/fs03/ag36/Shani/Trev-Seq/AnimalTFDB/Mus_musculus_TF.txt", header = TRUE, sep="\t")
# --------- DGE.15 --------------------

dge.15$genes$SYMBOL <- dge.15$genes$ENSEMBL
TF_15               <- TF_FULL[TF_FULL$Ensembl %in% dge.15$genes$SYMBOL,] 
TF_15               <- TF_FULL[TF_FULL$Ensembl %in% toptable.15[,1],] 
tfs                 <- c(TF_15$Ensembl)
# "Egr1" "Zim1"
rtni.15 <- tni.constructor(expData = cpm(dge.15, log=TRUE), 
                        regulatoryElements = tfs, 
                        rowAnnotation = dge.15$genes, 
                        colAnnotation = dge.15$samples)
# permutation
rtni.15 <- tni.permutation(rtni.15, nPermutations = 100)

# bootstrapping for confidence
rtni.15 <- tni.bootstrap(rtni.15)

# applies the ARACNe algorithm (Margolin, Nemenman, et al. 2006) to remove the weakest interaction in any triplet formed by two TFs
rtni.15 <- tni.dpi.filter(rtni.15)

tni.regulon.summary(rtni.15)

tni.regulon.summary(rtni.15, regulatoryElements = "Egr1")

regulons.15 <- tni.get(rtni.15, what = "regulons.and.mode", idkey = "SYMBOL")
head(regulons.15$EGR1)

g <- tni.graph(rtni.15, regulatoryElements = tfs)

pdf("./Regulons/reg_15.pdf")
    library(RedeR)
    rdp <- RedPort()
    calld(rdp)
    addGraph(rdp, g, layout= NULL)
    addLegend.color(rdp, g, type = "edge")
    addLegend.shape(rdp, g)
    relax(rdp, ps = TRUE)
dev.off()
