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
setwd("/fs03/ag36/Shani/Trev-Seq/featurecounts_edgeR_Mar_2023_dedup/P3")
here()
dir.create("./GeneSetTesting", showWarnings = FALSE)
source("/fs04/ag36/Shani/Trev-Seq/scripts/helper_scripts.R")

# --------------- Load data -----------------
dge.3    <- readRDS("./dge3.final.rds")
v.3      <- readRDS("./voom3.final.rds")
efit.3   <- readRDS("./efit3.final.rds")
design.3 <- readRDS("./design3.final.rds")

# ---------- Define contrast matrix ---------------
cont.matrix.3 <- makeContrasts(
  # at E0.5 
  E_L_3vsE_R_3 = E_L_3-E_R_3, 
  E_L_3vsC_3   = E_L_3-C_3, 
  E_R_3vsC_3   = E_R_3-C_3,
  E_L_3vsRest  = E_L_3-(E_R_3+C_3),
  E_R_3vsRest  = E_R_3-(E_L_3+C_3),
  E_3VsC_3     = (E_L_3+E_R_3)-C_3,
  levels=colnames(design.3)) 

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
dge.3 <- ensemble_to_entrez(dge.3)

# --------------------- Run mroast -----------------------------------------

# mroast on Hallmark geneset
mroast.MM(dge.3, efit.3, MM.h, design.3, "3", "hallmark", cont.matrix.3)

# mroast on Gene ontology
mroast.MM(dge.3, efit.3, MM.c5, design.3, "3", "GO", cont.matrix.3)

# mroast on KEGG
mroast.MM(dge.3, efit.3, MM.c2, design.3, "3", "KEGG", cont.matrix.3)

# mroast on Immunogenicity
mroast.MM(dge.3, efit.3, MM.c7, design.3, "3", "Immun", cont.matrix.3)
