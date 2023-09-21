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
#setwd("/fs03/ag36/Shani/Trev-Seq/featurecounts_edgeR_Mar_2023/E17")
setwd("/fs03/ag36/Shani/Trev-Seq/featurecounts_edgeR_Mar_2023_dedup/E17")
here()
dir.create("./GeneSetTesting", showWarnings = FALSE)
source("/fs04/ag36/Shani/Trev-Seq/scripts/helper_scripts.R")

# --------------- Load data -----------------
dge.17    <- readRDS("./dge17.final.rds")
v.17      <- readRDS("./voom17.final.rds")
efit.17   <- readRDS("./efit17.final.rds")
design.17 <- readRDS("./design17.final.rds")

# -------------- Define contrast matrix --------------
cont.matrix.17 <- makeContrasts(
  # at E17.5 
  E_L_17vsE_R_17 = E_L_17-E_R_17, 
  E_L_17vsC_17   = E_L_17-C_17, 
  E_R_17vsC_17   = E_R_17-C_17,
  E_L_17vsRest   = E_L_17-(E_R_17+C_17),
  E_R_17vsRest   = E_R_17-(E_L_17+C_17),
  E_17VsC_17     = (E_L_17+E_R_17)-C_17,
  levels=colnames(design.17)) 

# ---------------- Load metadata -----------------
# metadata <- readRDS("/fs03/ag36/Shani/Trev-Seq/featurecounts_edgeR/final_metadata_object.rds")
metadata <- readRDS("/fs03/ag36/Shani/Trev-Seq/featurecounts_edgeR_Mar_2023_dedup/metadata_dedup_full.rds")

# MSigDB files
MM.c2 <- readRDS("/fs03/ag36/Shani/MSigDB_v7.1/Mm.c2.all.v7.1.entrez.rds")
MM.c3 <- readRDS("/fs03/ag36/Shani/MSigDB_v7.1/Mm.c3.all.v7.1.entrez.rds")
MM.c4 <- readRDS("/fs03/ag36/Shani/MSigDB_v7.1/Mm.c4.all.v7.1.entrez.rds")
MM.c5 <- readRDS("/fs03/ag36/Shani/MSigDB_v7.1/Mm.c5.all.v7.1.entrez.rds")
MM.c6 <- readRDS("/fs03/ag36/Shani/MSigDB_v7.1/Mm.c6.all.v7.1.entrez.rds")
MM.c7 <- readRDS("/fs03/ag36/Shani/MSigDB_v7.1/Mm.c7.all.v7.1.entrez.rds")
MM.h <- readRDS("/fs03/ag36/Shani/MSigDB_v7.1/Mm.h.all.v7.1.entrez.rds")

# ------------------- Convert Ensemble to Entrez for DGE objects ----------
dge.17 <- ensemble_to_entrez(dge.17)

# --------------------- Run mroast -----------------------------------------

# mroast on Hallmark geneset
mroast.MM(dge.17, efit.17, MM.h, design.17, "17", "hallmark", cont.matrix.17)

# mroast on Gene ontology
mroast.MM(dge.17, efit.17, MM.c5, design.17, "17", "GO", cont.matrix.17)

# mroast on KEGG
mroast.MM(dge.17, efit.17, MM.c2, design.17, "17", "KEGG", cont.matrix.17)

# mroast on Immunogenicity
mroast.MM(dge.17, efit.17, MM.c7, design.17, "17", "Immun", cont.matrix.17)
