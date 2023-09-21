#' ----- Salmon-EdgeR pipeline on Trev-Seq dataset --------
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

# --------- Plot library size before normalisation ----------

tmp.df           <- as.data.frame(cbind(dge$samples$lib.size,rownames(dge$samples))) 
colnames(tmp.df) <- c("lib.size", "sample")
tmp.df$month     <- as.factor(substr(tmp.df$sample, 1,3))
  
pdf("./library.size.unnormalised.full.pdf")
  ggplot(tmp.df, aes(x    = as.factor(sample), 
                     y    = as.numeric(lib.size),
                     fill = month)) + 
    geom_bar(stat = "identity") +
    theme(legend.position="none") +
    scale_fill_manual(values = c("deepskyblue3","darkorange3","goldenrod1")) +
    facet_grid(~month, scales="free", space="free_x",) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ylab("library size")
dev.off()

pdf("./log.library.size.unnormalised.full.pdf")
ggplot(tmp.df, aes(x    = as.factor(sample), 
                   y    = log(as.numeric(lib.size)),
                   fill = month)) + 
  geom_bar(stat = "identity") +
  theme(legend.position="none") +
  scale_fill_manual(values = c("deepskyblue3","darkorange3","goldenrod1")) +
  facet_grid(~month, scales="free", space="free_x",) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("log(library size)")
dev.off()

# ----------- MDS plot before normalisation -------------

pdf("./mds.plot.unnormalised.expgroup.full.pdf", width = 10, height = 10)
plotMDS(dge,
        labels = dge$samples$samples.s,
        method ="logFC", 
        col = as.numeric(dge$samples$month), 
        cex = 0.5)
dev.off()

# pdf version by colouring based on the experimental Right and Left
pdf("./mds.plot.unnormalised.genotype_side.full.pdf", width = 10, height = 10)
plotMDS(dge,
        labels = dge$samples$samples.s,
        method ="logFC", 
        col = as.numeric(dge$samples$genotype.s), 
        cex = 0.5)
legend("bottomleft",
       legend = as.character(unique(dge$samples$genotype.s)),
       col = as.numeric(unique(dge$samples$genotype.s)),
       pch = 15,
       cex = 0.7,
       pt.cex = 0.5)
dev.off()

# Here there are four things to note
# 1. technical replicates: for September (2 lanes) tech reps need to be merged : edgeR::sumTechReps()
# 2. batch effect: batch effect is by the month. 
#           MDS plot shows this
#           We should ideally see all normalisation samples overlapping with eachother
#           But it is not the case here
#           So to remove the batch effect, we need to model the month variable in to the design matrix later
# 3. Remove the sheath samples from July from further analysis
#           These samples seem to be the of low quality (low coverage)
#           So they need to be taken out
# 4. Removing the normalisation samples from further analysis
#           They are so different to the rest, so even if we TMM normalise 
#           them, they will be not well normalised and will affect the normalisation factors

# getting a sense of the sizes in the dge object
L <- mean(dge$samples$lib.size) * 1e-6
M <- median(dge$samples$lib.size) * 1e-6
c(L, M)
# [1] 5.793607 4.979603

# ------------- Remove the Sheath samples ------------------
metadata <- metadata[metadata$BONE_LOC != "Sheath",]
dge      <- dge[,dge$samples$bone.loc != "Sheath"]

# ------------- Summing over the replicates -------------------
dge.colnames           <- colnames(dge) # to double check later
dge$samples$m.sample.s <- paste(dge$samples$month, dge$samples$sample.s, sep ="_")
dge                    <- sumTechReps(dge, ID=dge$samples$m.sample.s)

# adding a new column on metadata file for ordering
metadata$ID.NEW        <- factor(paste(metadata$MONTH, metadata$SAMPLE, sep ="_"))

# getting a sense of the sizes in the edited dge object
L <- mean(dge$samples$lib.size) * 1e-6
M <- median(dge$samples$lib.size) * 1e-6
c(L, M)
# [1] 8.053261 7.286095 # increased values

# --------------- Filter the DGElist object based on transcript counts ----------------
dim(dge)
# [1] 116758     56
dge.full             <- dge # keep the old one in case we mess up
head(dge$counts)
head(cpm(dge))

#apply(dge$counts, 2, sum) # total gene counts per sample

table(rowSums(dge$counts==0)==2)

keepMin <- apply(dge$counts, 1, max) >= 0.0
keepCpm <- rowSums(cpm(dge)> 0.5) >= 2                  # Keep only genes with cpm above x in at least y samples
keep    <- keepMin & keepCpm
dge     <- dge[keep,]

dim(dge)
# [1] 
dge$samples$lib.size <- colSums(dge$counts)


# getting a sense of the sizes in the edited dge object
L <- mean(dge$samples$lib.size) * 1e-6
M <- median(dge$samples$lib.size) * 1e-6
c(L, M)
# [1] 6.318618 5.476243

# ----------------- Create design matrix for batch effect detection ----------------------

# because we have merged the replicates,
# this has to reflect in the metadata as well
metadata.full  <- metadata
metadata       <- metadata[metadata$REPLICATE != 2, ]

group  <- as.factor(metadata$expgroup[match(colnames(dge$counts), metadata$ID.NEW)])
month  <- as.factor(metadata$MONTH[match(colnames(dge$counts), metadata$ID.NEW)])

design <- model.matrix(~0+group+month) 

# ----------------- Assesing the batch effect -----------------

tmp.dge  <- dge
tmp.dge  <- calcNormFactors(tmp.dge)
tmp.dge

dge.logcpm <- cpm(tmp.dge,log=TRUE,prior.count=2)
dge.logcpm <- removeBatchEffect(dge.logcpm, batch=as.character(tmp.dge$samples$month), design = design)

# MDS plot by month
pdf("./mds.plot.expgroup.batcheffectcorr.pdf", width = 10, height = 10)
    plotMDS(dge.logcpm, 
            labels = dge$samples$samples.s,
            col=as.numeric(dge$samples$month),
            cex = 0.5)
    legend("bottomleft",
        legend = as.character(unique(dge$samples$month)),
        col = as.numeric(unique(dge$samples$month)),
        pch = 15,
        cex = 0.7,
        pt.cex = 0.5)
dev.off()

# MDS plot by genotype and side
pdf("./mds.plot.genotype_side_batcheffectcorr.pdf", width = 10, height = 10)
    plotMDS(dge.logcpm,
            labels = dge$samples$samples.s,
            method ="logFC", 
            col = as.numeric(dge$samples$genotype.s), 
            cex = 0.5)
    legend("bottomleft",
        legend = as.character(unique(dge$samples$genotype.s)),
        col = as.numeric(unique(dge$samples$genotype.s)),
        pch = 15,
        cex = 0.7,
        pt.cex = 0.5)
dev.off()

# -------------- Demonstrate the issue with TMM normalisation w/ anchor samples in the dataset ------------

tmp.df           <- as.data.frame(cbind(tmp.dge$samples$norm.factors,rownames(tmp.dge$samples))) 
colnames(tmp.df) <- c("norm.factors", "sample")
tmp.df$month     <- as.factor(substr(tmp.df$sample, 1,3))
lcpm             <- cpm(tmp.dge, log=TRUE)

pdf("./norm.factors.normalised.withplacenta.pdf")
    ggplot(tmp.df, aes(x    = as.factor(sample), 
                    y    = as.numeric(norm.factors),
                    fill = month)) + 
    geom_bar(stat = "identity") +
    theme(legend.position="none") +
    scale_fill_manual(values = c("deepskyblue3","darkorange3","goldenrod1")) +
    facet_grid(~month, scales="free", space="free_x",) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ylab("normalisation factor values")
dev.off()

pdf("./norm.factors.boxplot.normalised.withplacenta.pdf")

   colourCount  <- length(unique(tmp.dge$samples$genotype.s))
   getPalette   <- colorRampPalette(brewer.pal(7, "Accent"))
   col_groups   <- getPalette(colourCount)

   par(cex.axis=0.45)
   boxplot(lcpm, las=2, col=col_groups[tmp.dge$samples$genotype.s], main="")
   title(ylab="Log-cpm")

dev.off()

rm(tmp.dge, tmp.df)

# ------------- Remove the Anchor (placenta) samples ------------------
metadata <- metadata[metadata$SAMPLE_NO != "Dahup",]
dge      <- dge[, dge$samples$sample != "Dahup"]

# getting a sense of the sizes in the dge object
L <- mean(dge$samples$lib.size) * 1e-6
M <- median(dge$samples$lib.size) * 1e-6
c(L, M)
# [1] 5.362832 5.381169 # with default filtering
# [1] 7.041437 6.953592 # with no filtering 

# -------------- TMM normalisation -----------------------

# The calcNormFactors() function normalizes for RNA composition by finding a set 
# of scaling factors for the library sizes that minimize the log-fold changes 
# between the samples for most genes. The default method for computing these 
# scale factors uses a trimmed mean of M-values (TMM) between each pair of 
# samples. We call the product of the original library size and the scaling 
# factor the effective library size.The effective library size replaces 
# the original library size in all downsteam analyses.
dge                  <- calcNormFactors(dge)
dge

# -------------- MDS plot for the normalised samples ------------
# pdf version for higher resolution
pdf("./mds.plot.normalised.expgroup.pdf", width = 10, height = 10)
plotMDS(dge,
        labels = dge$samples$samples.s,
        method ="logFC", 
        col = as.numeric(dge$samples$month), 
        cex = 0.9)
dev.off()

# pdf version by colouring based on the experimental Right and Left
pdf("./mds.plot.normalised.genotype_side.pdf", width = 10, height = 10)
   colourCount  <- length(unique(tmp.dge$samples$genotype.s))
   getPalette   <- colorRampPalette(brewer.pal(7, "Dark2"))
   col_groups   <- getPalette(colourCount)

plotMDS(dge,
        labels = dge$samples$samples.s,
        method ="logFC", 
        col = col_groups[dge$samples$genotype.s], 
        cex = 0.5)
legend("bottomleft",
       legend = as.character(unique(dge$samples$genotype.s)),
       col = col_groups[dge$samples$genotype.s],
       pch = 15,
       cex = 0.7,
       pt.cex = 0.5)
dev.off()

# ----------- Plotting the normalisatons factors for final dataset ----------

tmp.df           <- as.data.frame(cbind(dge$samples$norm.factors,rownames(dge$samples))) 
colnames(tmp.df) <- c("norm.factors", "sample")
tmp.df$month     <- as.factor(substr(tmp.df$sample, 1,3))
lcpm             <- cpm(dge, log=TRUE)

pdf("./norm.factors.normalised.nofiltering.pdf")
ggplot(tmp.df, aes(x    = as.factor(sample), 
                   y    = as.numeric(norm.factors),
                   fill = month)) + 
  geom_bar(stat = "identity") +
  theme(legend.position="none") +
  scale_fill_manual(values = c("deepskyblue3","darkorange3","goldenrod1")) +
  facet_grid(~month, scales="free", space="free_x",) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("normalisation factor values")
dev.off()

pdf("./norm.factors.boxplot.normalised.nofiltering.pdf")
   colourCount  <- length(unique(dge$samples$genotype.s))
   getPalette   <- colorRampPalette(brewer.pal(7, "Dark2"))
   col_groups   <- getPalette(colourCount)

   par(cex.axis=0.45)
   boxplot(lcpm, las=2, col=col_groups[dge$samples$genotype.s], main="")
   title(ylab="Log-cpm")
dev.off()

# ----------------- Density plot for the counts --------------------

# here tpm.dge is a normalised (TMM) dge object of 56 samples (excludiing the sheath)
lcpm <- cpm(dge,log=TRUE,prior.count=2)

pdf("./density.filtered.normalised.logcpm.no.filtering.pdf")

    lcpm.cutoff <- log2(10/M + 2/L)
    nsamples    <- ncol(dge)

    colourCount  <- length(unique(colnames(dge)))
    getPalette   <- colorRampPalette(brewer.pal(7, "Accent"))
    col_groups   <- getPalette(colourCount)

    plot(density(lcpm[,1]), col=col_groups[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
    title(xlab="Log-cpm")
    #abline(v=lcpm.cutoff, lty=3)
    for (i in 2:nsamples){
      den <- density(lcpm[,i])
      lines(den$x, den$y, col=col_groups[i], lwd=2)
    }
    legend("topright", 
    legend = as.character(unique(colnames(tmp.dge))),
    text.col=col_groups, bty="n",
    cex = 0.7)
dev.off()

# ----------------- Create design matrix for DE analysis ----------------------

group  <- as.factor(metadata$expgroup[match(colnames(dge$counts), metadata$ID.NEW)])
month  <- as.factor(metadata$MONTH[match(colnames(dge$counts), metadata$ID.NEW)])
#design <- model.matrix(~0+group+month) 
design <- model.matrix(~0+group)

# ----------------- Create BCV plot ------------------------------------------

dge <- estimateDisp(dge, design, robust=TRUE)
pdf("bcv.unfiltered.pdf")
    plotBCV(dge)
dev.off()

# -------------------- voom ---------------------

pdf("./voom.dge.normalised.pdf")
  v   <- voom(dge, design, plot = TRUE)
dev.off()

# --------------------- Fitting linear model --------------

#cor <- duplicateCorrelation(v, design, block = month)
#fit <- lmFit(object=v, design=design, 
#  correlation=cor$consensus.correlation)

fit <- lmFit(object=v, design=design)

# ------------------ contrast matrix --------------

cont.matrix <- makeContrasts(
  # at E15.5 
  E_L_15vsE_R_15 = groupE_L_15-groupE_R_15, 
  C_L_15vsC_R_15 = groupC_L_15-groupC_R_15,
  E_L_15vsC_L_15 = groupE_L_15-groupC_L_15, 
  E_R_15vsC_R_15 = groupE_R_15-groupC_R_15,
  # interaction
  E_15vsC_15     = (groupE_L_15-groupE_R_15)-(groupC_L_15 - groupC_R_15),

  # at 17.5
  E_L_17vsE_R_17 = groupE_L_17-groupE_R_17, 
  C_L_17vsC_R_17 = groupC_L_17-groupC_R_17,
  E_L_17vsC_L_17 = groupE_L_17-groupC_L_17, 
  E_R_17vsC_R_17 = groupE_R_17-groupC_R_17,
  # interaction
  E_17vsC_17     = (groupE_L_17-groupE_R_17)-(groupC_L_17 - groupC_R_17),
 
  # at P0
  E_L_0vsE_R_0 = groupE_L_0-groupE_R_0, 
  C_L_0vsC_R_0 = groupC_L_0-groupC_R_0,
  E_L_0vsC_L_0 = groupE_L_0-groupC_L_0, 
  E_R_0vsC_R_0 = groupE_R_0-groupC_R_0,
  # interaction
  E_0vsC_0     = (groupE_L_0-groupE_R_0)-(groupC_L_0 - groupC_R_0),

  # at P3
  E_L_3vsE_R_3 = groupE_L_3-groupE_R_3, 
  C_L_3vsC_R_3 = groupC_L_3-groupC_R_3,
  E_L_3vsC_L_3 = groupE_L_3-groupC_L_3, 
  E_R_3vsC_R_3 = groupE_R_3-groupC_R_3,
  # interaction
  E_3vsC_3     = (groupE_L_3-groupE_R_3)-(groupC_L_3 - groupC_R_3),
  levels=colnames(design)) 

# -------------------------- Contrast fitting --------------------

fit2 <- contrasts.fit(fit, cont.matrix)
efit <- eBayes(fit2, proportion=0.1)

summary(decideTests(efit), p.value=0.01,lfc=0)

#   E_L_15vsE_R_15 C_L_15vsC_R_15 E_L_15vsC_L_15 E_R_15vsC_R_15 E_15vsC_15
# Down              510             11              5             35        148
# NotSig         116155         116736         116729         116264     116595
# Up                 93             11             24            459         15
#        E_L_17vsE_R_17 C_L_17vsC_R_17 E_L_17vsC_L_17 E_R_17vsC_R_17 E_17vsC_17
# Down               39             99          69182             81         23
# NotSig         116662         116574          46395         116566     116704
# Up                 57             85           1181            111         31
#        E_L_0vsE_R_0 C_L_0vsC_R_0 E_L_0vsC_L_0 E_R_0vsC_R_0 E_0vsC_0
# Down             10            9           14            5        1
# NotSig       116741       116736       116712       116732   116755
# Up                7           13           32           21        2
#        E_L_3vsE_R_3 C_L_3vsC_R_3 E_L_3vsC_L_3 E_R_3vsC_R_3 E_3vsC_3
# Down             52        57319           14           13        5
# NotSig       116694        58868       116728       116720   116748
# Up               12          571           16           25        5


tfit <- treat(efit, lfc=1)
dt   <- decideTests(tfit)
summary(dt)

#           E_L_15vsE_R_15 C_L_15vsC_R_15 E_L_15vsC_L_15 E_R_15vsC_R_15 E_15vsC_15
# Down                0              0              0              0          0
# NotSig         116758         116758         116756         116758     116758
# Up                  0              0              2              0          0
#        E_L_17vsE_R_17 C_L_17vsC_R_17 E_L_17vsC_L_17 E_R_17vsC_R_17 E_17vsC_17
# Down                3              9              7             11          8
# NotSig         116749         116740         116737         116731     116734
# Up                  6              9             14             16         16
#        E_L_0vsE_R_0 C_L_0vsC_R_0 E_L_0vsC_L_0 E_R_0vsC_R_0 E_0vsC_0
# Down              0            5            0            1        0
# NotSig       116757       116752       116751       116755   116758
# Up                1            1            7            2        0
#        E_L_3vsE_R_3 C_L_3vsC_R_3 E_L_3vsC_L_3 E_R_3vsC_R_3 E_3vsC_3
# Down              0            0            0            0        0
# NotSig       116757       116758       116757       116756   116757
# Up                1            0            1            2        1

pdf("./efit.unfiltered.pdf")
  plotSA(efit)
dev.off()

out <- topTable(efit, n=Inf, adjust="BH")
out2 <- cbind(efit$coef,
              out[, c('P.Value','adj.P.Val','AveExpr')])

write.table(out2, "./toptable_all.csv", sep=",", quote = FALSE, row.names = TRUE)

glimma.res <- glMDPlot(efit, coef=1, status=dt, main=colnames(efit)[1],
         side.main="ENSEMBLID", counts=lcpm, groups=group, launch=FALSE)

htmlwidget::saveWidget(glimma.res, "./Glimma.total.unfiltered.html")

# add a look here to print all MD plots
pdf("mdplot_de.genes.pdf")
    #results <- as.data.frame(topTable(efit, n=Inf, adjust="BH"))
    plotMD(efit, column=1, status = decideTests(efit))
    #text(results$logCPM,results$logFC,labels = rownames(results),col="black",cex=0.5,pos=3)
dev.off()