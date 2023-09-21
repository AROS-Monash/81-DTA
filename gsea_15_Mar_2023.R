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
# module load R/4.2.2-mkl

# ------------- Set working directory ----
# location of processing
setwd("/fs03/ag36/Shani/Trev-Seq/featurecounts_edgeR_Mar_2023_dedup/E15")
here::here()
dir.create("./GeneSetTesting", showWarnings = FALSE)
source("/fs04/ag36/Shani/Trev-Seq/scripts/helper_scripts.R")

# --------------- Load data -----------------
dge.15    <- readRDS("./dge15.final.rds")
v.15      <- readRDS("./voom15.final.rds")
efit.15   <- readRDS("./efit15.final.rds")
design.15 <- readRDS("./design15.final.rds")

# -------------- Define contrast matrix --------------
cont.matrix.15 <- makeContrasts(
  # at E15.5 
  E_L_15vsE_R_15 = E_L_15-E_R_15, 
  E_L_15vsC_15   = E_L_15-C_15, 
  E_R_15vsC_15   = E_R_15-C_15,
  E_L_15vsRest   = E_L_15-(E_R_15+C_15),
  E_R_15vsRest   = E_R_15-(E_L_15+C_15),
  E_15VsC_15     = (E_L_15+E_R_15)-C_15,
  #EVsC           = (E_L_15-C_15)-(E_R_15-C_15), # same as E_L vs E_R
  levels=colnames(design.15)) 

# ---------------- Load metadata -----------------
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
dge.15 <- ensemble_to_entrez(dge.15)

# --------------------- Run mroast -----------------------------------------

# mroast on Hallmark geneset
mroast.MM(dge.15, efit.15, MM.h, design.15, "15", "hallmark", cont.matrix.15)

# mroast on Gene ontology
mroast.MM(dge.15, efit.15, MM.c5, design.15, "15", "GO", cont.matrix.15)

# mroast on KEGG
mroast.MM(dge.15, efit.15, MM.c2, design.15, "15", "KEGG", cont.matrix.15)

# mroast on Immunogenicity
mroast.MM(dge.15, efit.15, MM.c7, design.15, "15", "Immun", cont.matrix.15)
