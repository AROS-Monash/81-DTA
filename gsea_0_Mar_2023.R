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
setwd("/fs03/ag36/Shani/Trev-Seq/featurecounts_edgeR_Mar_2023_dedup/P0")
here()
dir.create("./GeneSetTesting", showWarnings = FALSE)
source("/fs04/ag36/Shani/Trev-Seq/scripts/helper_scripts.R")

# --------------- Load data -----------------
dge.0    <- readRDS("./dge0.final.rds")
v.0      <- readRDS("./voom0.final.rds")
efit.0   <- readRDS("./efit0.final.rds")
design.0 <- readRDS("./design0.final.rds")

# -------------- Define contrast matrix --------------
cont.matrix.0 <- makeContrasts(
  # at E0.5 
  E_L_0vsE_R_0 = E_L_0-E_R_0, 
  E_L_0vsC_0   = E_L_0-C_0, 
  E_R_0vsC_0   = E_R_0-C_0,
  E_L_0vsRest  = E_L_0-(E_R_0+C_0),
  E_R_0vsRest  = E_R_0-(E_L_0+C_0),
  E_0VsC_0     = (E_L_0+E_R_0)-C_0,
  levels=colnames(design.0)) 

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
dge.0 <- ensemble_to_entrez(dge.0)

# --------------------- Run mroast -----------------------------------------

# mroast on Hallmark geneset
mroast.MM(dge.0, efit.0, MM.h, design.0, "0", "hallmark", cont.matrix.0)

# mroast on Gene ontology
mroast.MM(dge.0, efit.0, MM.c5, design.0, "0", "GO", cont.matrix.0)

# mroast on KEGG
mroast.MM(dge.0, efit.0, MM.c2, design.0, "0", "KEGG", cont.matrix.0)

# mroast on Immunogenicity
mroast.MM(dge.0, efit.0, MM.c7, design.0, "0", "Immun", cont.matrix.0)
