#' ----- EdgeR-Limma pipeline on Trev-Seq dataset (transcript level) --------
#' @author Shani Amaraisnghe
#' @date 11/07/2022
#' @aliases Rosello-Diez lab, ARMI, Monash
#' @concept analysing Trev-seq bulk RNA-Seq data for Chee
#' @description for the july-2021 the data was demultiplexed using umi-tools, and bbmap 
#' had to spend two rounds of demultiplexing to get the data seperated to samples and the UMIs.
#' For sept-2021 data the original fastq files available in the raw data directory was used as 
#' they were already demultiplexed
#' for nov-2021 data, one round of umi-tools was carried out to demultiplex by UMI  
#' next, the demux data was aligned and quantified using Salmon
#' now, need to look at the them via edgeR
#' running the edgeR-limma pipeline in M3
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

# -------------I/O --------------
# working in M3
#module load R/4.2-mkl

# ------------- Set working directory ----

# location of processing
setwd("/fs03/ag36/Shani/Trev-Seq/salmon_edgeR")
here()

# ------------- Load the DGEList object ------------

# the script to create the DGE list object is in the salmon_to_dge.R
dge <- readRDS("/fs03/ag36/Shani/Trev-Seq/salmon-edgeR/complete_DGEObject_unnormalised.RDS")

# ------------ Load metadata ---------------

metadata <- read.table("/fs02/ag36/Shani/Trev-Seq/fastq_exp.csv", sep = ",", header = TRUE)

# --------- Modifying metadata -------------

metadata$expgroup<- paste(substr(metadata$GENOTYPE, 1,1),
                            metadata$SIDE,
                            substr(metadata$STAGE, 2,3),
                            sep = "_")

# ------------- Modifying dge$samples -------------

dge$samples$sample        <- metadata$SAMPLE_NO[match(row.names(dge$samples),metadata$ID)]
dge$samples$month         <- as.factor(metadata$MONTH[match(row.names(dge$samples),metadata$ID)])
dge$samples$stage         <- as.factor(metadata$STAGE[match(row.names(dge$samples),metadata$ID)])
dge$samples$genotype      <- as.factor(metadata$GENOTYPE[match(row.names(dge$samples),metadata$ID)])
dge$samples$side          <- as.factor(metadata$SIDE[match(row.names(dge$samples),metadata$ID)])
dge$samples$replicate     <- metadata$REPLICATE[match(row.names(dge$samples),metadata$ID)]
dge$samples$expgroup      <- as.factor(metadata$expgroup[match(row.names(dge$samples),metadata$ID)])
dge$samples$sample.s      <- as.factor(metadata$SAMPLE[match(row.names(dge$samples),metadata$ID)])
dge$samples$genotype.s    <- as.factor(paste(dge$samples$genotype, dge$samples$side, sep = "_"))
dge$samples$bone.loc      <- as.factor(metadata$BONE_LOC[match(row.names(dge$samples),metadata$ID)])

# ------------- Remove the Sheath samples ------------------
metadata <- metadata[metadata$BONE_LOC != "Sheath",]
dge      <- dge[,dge$samples$bone.loc != "Sheath"]

# ------------- Summing over the replicates -------------------
dge.colnames           <- colnames(dge) # to double check later
dge$samples$m.sample.s <- paste(dge$samples$month, dge$samples$sample.s, sep ="_")
dge                    <- sumTechReps(dge, ID=dge$samples$m.sample.s)

# ------------- Remove the Anchor (placenta) samples ------------------
metadata <- metadata[metadata$SAMPLE_NO != "Dahup",]
dge      <- dge[, dge$samples$sample != "Dahup"]

# adding a new column on metadata file for ordering
metadata$ID.NEW        <- factor(paste(metadata$MONTH, metadata$SAMPLE, sep ="_"))

############################################################
# --------------------- Stage: E 15.5 -----------------------
############################################################

dge.15 <- dge[,as.character(metadata$ID.NEW[metadata$STAGE == "E15_5"])]

# filtering

dge.15.full  <- dge.15 # keep the old one in case we mess up
head(dge.15$counts)
head(cpm(dge.15))

#apply(dge.15$counts, 1, max) # total gene counts per sample

table(rowSums(dge.15$counts==0)==2)

keepMin <- apply(dge.15$counts, 1, max) >= 0.0
keepCpm <- rowSums(cpm(dge.15)>1) >= 2                  # Keep only genes with cpm above x in at least y samples
keep    <- keepMin & keepCpm
dge.15  <- dge.15[keep,]

# keep    <- rowSums(cpm(dge.15)>100) >= 2
# dge.15  <- dge.15[keep,]

dim(dge.15)
# [1] 31754    20

dge.15$samples$lib.size <- colSums(dge.15$counts)

# TMM normalisation

dge.15   <- calcNormFactors(dge.15)
dge.15

# design matrix

group.15            <- as.factor(metadata$expgroup[match(colnames(dge.15$counts), metadata$ID.NEW)])
month.15            <- as.factor(metadata$MONTH[match(colnames(dge.15$counts), metadata$ID.NEW)])
design.15           <- model.matrix(~0+group.15+month.15)
colnames(design.15) <- gsub("group\\.15","", colnames(design.15))

# voom

pdf("./voom.dge.15.normalised.pdf")
  v.15   <- voom(dge.15, design.15, plot = TRUE)
dev.off()

# Fitting linear model

#cor.15 <- duplicateCorrelation(v.15, design.15, block = month.15)
fit.15 <- lmFit(object=v.15, design=design.15)

#  contrast matrix 

cont.matrix.15 <- makeContrasts(
  # at E15.5 
  E_L_15vsE_R_15 = E_L_15-E_R_15, 
  C_L_15vsC_R_15 = C_L_15-C_R_15,
  E_L_15vsC_L_15 = E_L_15-C_L_15, 
  E_R_15vsC_R_15 = E_R_15-C_R_15,
  # interaction
  E_15vsC_15     = (E_L_15-E_R_15)-(C_L_15 - C_R_15),
  levels=colnames(design.15)) 

# contrast fitting

fit2.15 <- contrasts.fit(fit.15, cont.matrix.15)
efit.15 <- eBayes(fit2.15)

summary(decideTests(efit.15), p.value=0.01,lfc=0)

# E_L_15vsE_R_15 C_L_15vsC_R_15 E_L_15vsC_L_15 E_R_15vsC_R_15 E_15vsC_15
# Down                0              0              0              0          0
# NotSig          31754          31754          31750          31752      31754
# Up                  0              0              4              2          0

pdf("./efit.15.pdf")
  plotSA(efit.15)
dev.off()

# DE results

out.15 <- topTable(efit.15, n=Inf, sort.by='F', adjust="BH")
out2.15 <- cbind(efit.15$coef,
              out.15[, c('P.Value','adj.P.Val','AveExpr')])

write.table(out2.15, "./toptable_15.csv", sep=",", quote = FALSE, row.names = TRUE)

############################################################
# --------------------- Stage: E 17.5 -----------------------
############################################################

dge.17 <- dge[,as.character(metadata$ID.NEW[metadata$STAGE == "E17_5"])]

# filtering 

dge.17.full  <- dge.17 # keep the old one in case we mess up
head(dge.17$counts)
head(cpm(dge.17))

#apply(dge.17$counts, 1, max) # total gene counts per sample

table(rowSums(dge.17$counts==0)==2)

keepMin <- apply(dge.17$counts, 1, max) >= 0.0
keepCpm <- rowSums(cpm(dge.17)>0) >= 3                  # Keep only genes with cpm above x in at least y samples
keep    <- keepMin & keepCpm
dge.17  <- dge.17[keep,]

dim(dge.17)
# [1]  32518    10
dge.17$samples$lib.size <- colSums(dge.17$counts)

# TMM normalisation

dge.17 <- calcNormFactors(dge.17)
dge.17

# design matrix

group.17            <- as.factor(metadata$expgroup[match(colnames(dge.17$counts), metadata$ID.NEW)])
design.17           <- model.matrix(~0+group.17)
colnames(design.17) <- gsub("group\\.17","", colnames(design.17))
design.17

# voom

pdf("./voom.dge.17.normalised.pdf")
  v.17   <- voom(dge.17, design.17, plot = TRUE)
dev.off()

# fitting linear model

#cor.17 <- duplicateCorrelation(v.17, design.17)
#fit.17 <- lmFit(object=v.17, design=design.17, 
#  correlation=cor.17$consensus.correlation)

fit.17 <- lmFit(object=v.17, design=design.17)

#  contrast matrix 

cont.matrix.17 <- makeContrasts(
  # at 17.5
  E_L_17vsE_R_17 = E_L_17-E_R_17, 
  C_L_17vsC_R_17 = C_L_17-C_R_17,
  E_L_17vsC_L_17 = E_L_17-C_L_17, 
  E_R_17vsC_R_17 = E_R_17-C_R_17,
  # interaction
  E_17vsC_17     = (E_L_17-E_R_17)-(C_L_17 - C_R_17),
  levels=colnames(design.17)) 

#  contrast fitting 

fit2.17 <- contrasts.fit(fit.17, cont.matrix.17)
efit.17 <- eBayes(fit2.17)

summary(decideTests(efit.17))

# no filtering
#        E_L_17vsE_R_17 C_L_17vsC_R_17 E_L_17vsC_L_17 E_R_17vsC_R_17 E_17vsC_17
# Down               25             38             42             69         12
# NotSig         116686         116687         116663         116630     116735
# Up                 47             33             53             59         11

# filtering
#        E_L_17vsE_R_17 C_L_17vsC_R_17 E_L_17vsC_L_17 E_R_17vsC_R_17 E_17vsC_17
# Down                9             10             11             33          0
# NotSig          37281          37292          37274          37240      37307
# Up                 17              5             22             34          0

#newfiltering
#        E_L_17vsE_R_17 C_L_17vsC_R_17 E_L_17vsC_L_17 E_R_17vsC_R_17 E_17vsC_17
# Down                0              0              0              0          0
# NotSig          26913          26913          26913          26913      26913
# Up                  0              0              0              0          0

# current filtering
#        E_L_17vsE_R_17 C_L_17vsC_R_17 E_L_17vsC_L_17 E_R_17vsC_R_17 E_17vsC_17
# Down                3              0              4              4          0
# NotSig          32506          32518          32512          32509      32518
# Up                  9              0              2              5          0

pdf("./efit.17.pdf")
  plotSA(efit.17)
dev.off()

# DE results

out.17 <- topTable(efit.17, n=Inf, sort.by='F', adjust="BH")
out2.17 <- cbind(efit.17$coef,
              out.17[, c('P.Value','adj.P.Val','AveExpr')])

# the unfiltered results are saved
write.table(out2.17, "./toptable_17.csv", sep=",", quote = FALSE, row.names = TRUE)

############################################################
# ---------------------   Stage: P0  -----------------------
############################################################

dge.0 <- dge[,as.character(metadata$ID.NEW[metadata$STAGE == "P0"])]

# filtering

dge.0.full  <- dge.0 # keep the old one in case we mess up
head(dge.0$counts)
head(cpm(dge.0))

table(rowSums(dge.0$counts==0)==2)

keepMin <- apply(dge.0$counts, 1, max) >= 0.0
keepCpm <- rowSums(cpm(dge.0)> 1) >= 10                  # Keep only genes with cpm above x in at least y samples
keep    <- keepMin & keepCpm
dge.0  <- dge.0[keep,]

dim(dge.0)
# [1] 21492    20
dge.0$samples$lib.size <- colSums(dge.0$counts)

# TMM normalisation

dge.0 <- calcNormFactors(dge.0)
dge.0

# design matrix

group.0            <- as.factor(metadata$expgroup[match(colnames(dge.0$counts), metadata$ID.NEW)])
design.0           <- model.matrix(~0+group.0)
colnames(design.0) <- gsub("group\\.0","", colnames(design.0))
design.0

# voom

pdf("./voom.dge.0.normalised.pdf")
  v.0   <- voom(dge.0, design.0, plot = TRUE)
dev.off()

# fitting linear model

fit.0 <- lmFit(object=v.0, design=design.0)

#  contrast matrix 

cont.matrix.0 <- makeContrasts(
  # at P0
  E_L_0vsE_R_0 = E_L_0-E_R_0, 
  C_L_0vsC_R_0 = C_L_0-C_R_0,
  E_L_0vsC_L_0 = E_L_0-C_L_0, 
  E_R_0vsC_R_0 = E_R_0-C_R_0,
  # interaction
  E_0vsC_0     = (E_L_0-E_R_0)-(C_L_0 - C_R_0),
  levels=colnames(design.0)) 

#  contrast fitting 

fit2.0 <- contrasts.fit(fit.0, cont.matrix.0)
efit.0 <- eBayes(fit2.0)

summary(decideTests(efit.0), p.value=0.01,lfc=0)

#        E_L_0vsE_R_0 C_L_0vsC_R_0 E_L_0vsC_L_0 E_R_0vsC_R_0 E_0vsC_0
# Down            661          449         1077          908      412
# NotSig        20301        20501        19485        19762    20752
# Up              530          542          930          822      328

pdf("./efit.0.pdf")
  plotSA(efit.0)
dev.off()

# DE results

out.0 <- topTable(efit.0, n=Inf, sort.by='F', adjust="BH")
out2.0 <- cbind(efit.0$coef,
              out.0[, c('P.Value','adj.P.Val','AveExpr')])

write.table(out2.0, "./toptable_0.csv", sep=",", quote = FALSE, row.names = TRUE)


############################################################
# ---------------------   Stage: P3  -----------------------
############################################################

dge.3 <- dge[,as.character(metadata$ID.NEW[metadata$STAGE == "P3"])]

# filtering

dge.3.full  <- dge.3 # keep the old one in case we mess up
head(dge.3$counts)
head(cpm(dge.3))

table(rowSums(dge.3$counts==3)==2)

keepMin <- apply(dge.3$counts, 1, max) >= 0.0
keepCpm <- rowSums(cpm(dge.3)> 1) >= 10                  # Keep only genes with cpm above x in at least y samples
keep    <- keepMin & keepCpm
dge.3  <- dge.3[keep,]

dim(dge.3)
# [1] 20554    20
dge.3$samples$lib.size <- colSums(dge.3$counts)

# TMM normalisation

dge.3 <- calcNormFactors(dge.3)
dge.3

# design matrix

group.3            <- as.factor(metadata$expgroup[match(colnames(dge.3$counts), metadata$ID.NEW)])
design.3           <- model.matrix(~0+group.3)
colnames(design.3) <- gsub("group\\.3","", colnames(design.3))
design.3

# voom

pdf("./voom.dge.3.normalised.pdf")
  v.3   <- voom(dge.3, design.3, plot = TRUE)
dev.off()

# fitting linear model

fit.3 <- lmFit(object=v.3, design=design.3)

#  contrast matrix 

cont.matrix.3 <- makeContrasts(
  # at P3
  E_L_3vsE_R_3 = E_L_3-E_R_3, 
  C_L_3vsC_R_3 = C_L_3-C_R_3,
  E_L_3vsC_L_3 = E_L_3-C_L_3, 
  E_R_3vsC_R_3 = E_R_3-C_R_3,
  # interaction
  E_3vsC_3     = (E_L_3-E_R_3)-(C_L_3 - C_R_3),
  levels=colnames(design.3)) 

#  contrast fitting 

fit2.3 <- contrasts.fit(fit.3, cont.matrix.3)
efit.3 <- eBayes(fit2.3)

summary(decideTests(efit.3), p.value=0.01,lfc=0)

#        E_L_3vsE_R_3 C_L_3vsC_R_3 E_L_3vsC_L_3 E_R_3vsC_R_3 E_3vsC_3
# Down            468          671          462          770      439
# NotSig        19405        19290        19629        19270    19579
# Up              681          593          463          514      536

pdf("./efit.3.pdf")
  plotSA(efit.3)
dev.off()

# DE results

out.3 <- topTable(efit.3, n=Inf, sort.by='F', adjust="BH")
out2.3 <- cbind(efit.3$coef,
              out.3[, c('P.Value','adj.P.Val','AveExpr')])

write.table(out2.3, "./toptable_3.csv", sep=",", quote = FALSE, row.names = FALSE)

lcpm.3  <- cpm(dge.3, log=TRUE)
dt.3    <- decideTests(efit.3)
glMDPlot(efit.3, coef=1, status=dt.3, main=colnames(efit.3)[1],
         side.main="ENSEMBLID", counts=lcpm.3, groups=group.3, launch=FALSE)

pdf("mdplot_de.genes.3.pdf", onefile=TRUE)
    for (i in length(colnames(efit.3))) {
        #results <- as.data.frame(topTable(efit, n=Inf, adjust="BH"))
        print(plotMD(efit.3, column=i, status = decideTests(efit.3)))
        #text(results$logCPM,results$logFC,labels = rownames(results),col="black",cex=0.5,pos=3)
    }
dev.off()
