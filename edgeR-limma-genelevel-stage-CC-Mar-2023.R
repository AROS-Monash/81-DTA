#' @author Shani Amarsinghe
#' @date 01/02/2023
#' @aliases Rosello-Diez lab, ARMI, Monash
#' @concept analysing Trev-seq bulk RNA-Seq data for Chee
#' @description for the july-2021 the data was demultiplexed using umi-tools, and bbmap 
#' had to spend two rounds of demultiplexing to get the data seperated to samples and the UMIs.
#' For sept-2021 data the original fastq files available in the raw data directory was used as 
#' they were already demultiplexed
#' for nov-2021 data, one round of umi-tools was carried out to demultiplex by UMI  
#' next, the demux data was aligned and read counting was done using featureCounts()
#' now, need to look at the them via edgeR
#' running the edgeR-limma pipeline in M3 for each stage (no July)
#' combining the Control samples
#' removing the some controls for E15.5
#' THIS IS DONE IN MARCH 2023 FOR THE LESS STRICT ALIGNMENT (MORE MISSMATCHES AND MORE MULTIMAPPING)
#' PURPOSE: TO UNCOVER TNA ALPHA IN GSEA AGAIN
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
library(plyr)
library(Glimma)
library(RColorBrewer)
library(Mus.musculus)
library(EnhancedVolcano, lib.loc="/fs03/ag36/Shani/R_libs/")
library(pheatmap)

# -------------I/O --------------
# working in M3
# module load R/4.2.2-mkl

# ------------- Set working directory ----
dir.create("/fs03/ag36/Shani/Trev-Seq/featurecounts_edgeR_Mar_2023", showWarnings = FALSE)

setwd("/fs03/ag36/Shani/Trev-Seq/featurecounts_edgeR_Mar_2023")
here::here()
source("/fs04/ag36/Shani/Trev-Seq/scripts/helper_scripts.R")

# ------------ I/O Read in feature counts (MM) ------------
counts.mm   <- read.table("/fs03/ag36/Shani/Trev-Seq/Counts_matrix/trev-seq-mm-countmatrix_Mar2023.tab")


# ------------- I/O Read in the metadata file -----------------
metadata <- read.table("/fs04/ag36/Shani/Trev-Seq/fastq_exp.csv", sep = ",", header = TRUE)

metadata$expgroup <- paste(substr(metadata$GENOTYPE, 1,1),
                            metadata$SIDE,
                            substr(metadata$STAGE, 2,3),
                            sep = "_")

metadata$ID.NEW  <- factor(paste(metadata$MONTH, metadata$SAMPLE, sep ="_"))

metadata$ID.2    <- factor(paste(metadata$MONTH, metadata$SAMPLE, metadata$REPLICATE, sep ="_"))

# ------------- Remove July samples from metadata ------
metadata <- metadata[metadata$MONTH != "Jul",]

# -------------- Adjust the column names for creating the DGElist object -----------
colnames(counts.mm)   <- metadata$ID.2[match(colnames(counts.mm), metadata$BAM_LOC)]
dge                   <- DGEList(counts = counts.mm, group = as.factor(metadata$expgroup[match(colnames(counts.mm), metadata$ID.2)]))

# ------------- Modifying dge$samples -------------
dge$samples$sample        <- metadata$SAMPLE_NO[match(row.names(dge$samples),metadata$ID.2)]
dge$samples$month         <- as.factor(metadata$MONTH[match(row.names(dge$samples),metadata$ID.2)])
dge$samples$stage         <- as.factor(metadata$STAGE[match(row.names(dge$samples),metadata$ID.2)])
dge$samples$genotype      <- as.factor(metadata$GENOTYPE[match(row.names(dge$samples),metadata$ID.2)])
dge$samples$side          <- as.factor(metadata$SIDE[match(row.names(dge$samples),metadata$ID.2)])
dge$samples$replicate     <- metadata$REPLICATE[match(row.names(dge$samples),metadata$ID.2)]
dge$samples$expgroup      <- as.factor(metadata$expgroup[match(row.names(dge$samples),metadata$ID.2)])
dge$samples$sample.s      <- as.factor(metadata$SAMPLE[match(row.names(dge$samples),metadata$ID.2)])
dge$samples$genotype.s    <- as.factor(paste(dge$samples$genotype, dge$samples$side, sep = "_"))
dge$samples$bone.loc      <- as.factor(metadata$BONE_LOC[match(row.names(dge$samples),metadata$ID.2)])

# ------------- Remove the Sheath samples ------------------
metadata <- metadata[metadata$BONE_LOC != "Sheath",]
dge      <- dge[,dge$samples$bone.loc != "Sheath"]

# ------------- Summing over the replicates -------------------
dge.colnames           <- colnames(dge) # to double check later
dge$samples$m.sample.s <- paste(dge$samples$month, dge$samples$sample.s, sep ="_")
dge                    <- sumTechReps(dge, ID=dge$samples$m.sample.s)

# because we have merged the replicates,
# this has to reflect in the metadata as well
metadata.full  <- metadata
metadata       <- metadata[metadata$REPLICATE != 2, ]

# ------------- Remove the Anchor (placenta) samples ------------------
metadata <- metadata[metadata$SAMPLE_NO != "Dahup",]
dge      <- dge[, dge$samples$sample != "Dahup"]

# ------- need to add an extra column to metadata to combine the controls -------
metadata$expgroup <- ifelse(metadata$GENOTYPE == "Ctrl",
                            paste(substr(metadata$GENOTYPE, 1,1),
                            substr(metadata$STAGE, 2,3),
                            sep = "_"),
                            paste(substr(metadata$GENOTYPE, 1,1),
                            metadata$SIDE,
                            substr(metadata$STAGE, 2,3),
                            sep = "_")
)

dge$samples$expgroup       <- as.factor(metadata$expgroup[match(row.names(dge$samples),metadata$ID.NEW)])

# -------------- Extracting common gene names from Ensembl for annotations of DE genes --------
genes <- select(Mus.musculus, keys=row.names(dge), columns=c("SYMBOL", "TXCHROM", "ENSEMBL"), 
                keytype="ENSEMBL")
dim(genes)

dge$genes <- as.data.frame(
    cbind(ENSEMBL = row.names(dge), 
    COMMON = ifelse(
        !is.na(
            genes$SYMBOL[match(row.names(dge), genes$ENSEMBL)]),
                            genes$SYMBOL[match(row.names(dge), genes$ENSEMBL)],
                            row.names(dge$counts)))
)
rownames(dge$genes) <- dge$genes$ENSEMBL

saveRDS(dge, "./dge.full.final.rds")

############################################################
# --------------------- Stage: E 15.5 -----------------------
############################################################

dir.create("/fs03/ag36/Shani/Trev-Seq/featurecounts_edgeR_Mar_2023/E15", showWarnings = FALSE)
setwd("/fs03/ag36/Shani/Trev-Seq/featurecounts_edgeR_Mar_2023/E15")

dge.15 <- dge[,as.character(metadata$ID.NEW[metadata$STAGE == "E15_5" & metadata$MONTH == "Nov" & metadata$SAMPLE_NO != "677"])]

# filtering
dge.15.full  <- dge.15 # keep the old one in case we mess up
head(dge.15$counts)
head(cpm(dge.15))

table(rowSums(dge.15$counts==0)==2)

# filter 2
group  <- as.factor(metadata$expgroup[match(colnames(dge.15$counts), metadata$ID.NEW)])

keep   <- filterByExpr(dge.15, group=group)
dge.15 <- dge.15[keep,]

dim(dge.15)
# [1] 20341    10 mm, filter 2, 2023 alignment, no july, no sample 677, combined controls

dge.15$samples$lib.size <- colSums(dge.15$counts)

# TMM normalisation
dge.15   <- calcNormFactors(dge.15)
dge.15

# normalisation factors
tmp.df           <- as.data.frame(cbind(dge.15$samples$norm.factors,rownames(dge.15$samples))) 
colnames(tmp.df) <- c("norm.factors", "sample")
tmp.df$month     <- as.factor(substr(tmp.df$sample, 1,3))
lcpm             <- cpm(dge.15, log=TRUE)

pdf("./norm.factors.normalised.15.pdf")
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

pdf("./norm.factors.boxplot.normalised.15.pdf")
   colourCount  <- length(unique(dge.15$samples$genotype.s))
   getPalette   <- colorRampPalette(brewer.pal(7, "Dark2"))
   col_groups   <- getPalette(colourCount)

   par(cex.axis=0.45)
   boxplot(lcpm, las=2, col=col_groups[dge.15$samples$genotype.s], main="")
   title(ylab="Log-cpm")
dev.off()

# MDS plot
pdf("./mds.plot.15.pdf", width = 10, height = 10)
    plotMDS(dge.15, 
            labels = dge.15$samples$m.samples.s,
            col=as.numeric(dge.15$samples$genotype),
            cex = 0.5)
    legend("bottomleft",
        legend = as.character(unique(dge.15$samples$genotype)),
        col = as.numeric(unique(dge.15$samples$genotype)),
        pch = 15,
        cex = 0.7,
        pt.cex = 0.5)
dev.off()

# design matrix
group.15            <- as.factor(metadata$expgroup[match(colnames(dge.15$counts), metadata$ID.NEW)])
design.15           <- model.matrix(~0+group.15)
colnames(design.15) <- gsub("group\\.15","", colnames(design.15))

# adding surrogate variables
sv <- wsva(cpm(dge.15, log=TRUE, prior.count=2),
        design=design.15, n.sv = 1, 
        weight.by.sd = TRUE, plot = FALSE)

design.15           <- model.matrix(~0 + group.15 + sv)
colnames(design.15) <- gsub("group\\.15","", colnames(design.15))

# BCV plot
dge.15 <- estimateDisp(dge.15, design.15, robust=TRUE)
pdf("bcv.15.pdf")
    plotBCV(dge.15)
dev.off()

# MDS plot before voom fitting
pdf("./mds.plot.15.final.pdf", width = 10, height = 10)
    plotMDS(dge.15, 
            labels = dge.15$samples$m.samples.s,
            col=as.numeric(dge.15$samples$genotype),
            cex = 0.5)
    legend("bottomleft",
        legend = as.character(unique(dge.15$samples$genotype)),
        col = as.numeric(unique(dge.15$samples$genotype)),
        pch = 15,
        cex = 0.7,
        pt.cex = 0.5)
dev.off()

# voom
pdf("./voom.dge.15.normalised.pdf")
  v.15   <- voomWithQualityWeights(dge.15, design.15, plot = TRUE)
dev.off()

# Fitting linear model
fit.15 <- lmFit(object=v.15, design=design.15)

#  contrast matrix 
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

# contrast fitting
fit2.15           <- contrasts.fit(fit.15, cont.matrix.15)
efit.15           <- eBayes(fit2.15)

# look at DE results
results.15 <- decideTests(efit.15, 
                       method = "separate", 
                       adjust.method = "fdr",
                       p.value=0.05,
                       lfc=0)

print(summary(results.15))

#        E_L_15vsE_R_15 E_L_15vsC_15 E_R_15vsC_15 E_L_15vsRest E_R_15vsRest
# Down                0           10            1        13445        13579
# NotSig          20340        20317        20335         6817         6735
# Up                  1           14            5           79           27
#        E_15VsC_15
# Down           29
# NotSig       7019
# Up          13293

write.table(summary(results.15),
            file      = "./de_summary_15_p1.txt",
            quote     = FALSE,
            sep       = "\t",
            row.names = TRUE)

pdf("./efit.15.pdf")
  plotSA(efit.15)
dev.off()

# save the v.15, dge.15 and efit.15
saveRDS(dge.15, "./dge15.final.rds")
saveRDS(v.15, "./voom15.final.rds")
saveRDS(efit.15, "./efit15.final.rds")
saveRDS(design.15, "./design15.final.rds")

# top tables

out.15           <- topTable(efit.15, adjust="BH", p.value=0.1, lfc=log2(1), n = Inf, sort.by = "F")
write.table(out.15, "./toptable_15_p1.csv", sep=",", quote = FALSE, row.names = TRUE)

# seeing a lot of outliers in the BCV plot - map RPKM values to find outliers and culprits of the high BCV, not attepted yet.

# MD plot

pdf("./MD_plots_all_p1_15.pdf")
    for (coef in colnames(efit.15)) {
        plotMD(efit.15, 
                coef = coef, 
                main = coef, 
                cex = c(0.3,0.5,0.5), 
                tag.genes = rownames(efit.15$genes)%in% row.names(results.15[which(results.15[,coef]!=0),])
        )
    }
dev.off()

# interactive plots - MD

for (coef.id in 1:length(colnames(efit.15))) {
    md_plot_interactive(efit.15, dge.15, coef.id, colnames(efit.15)[coef.id]) 
}

# volcano plot

pdf("./volcano_plots_all_p1_15.pdf")
    for (coef in colnames(efit.15)) {
        checkpoint <- row.names(results.15[which(results.15[,coef]!=0),])
        if (length(checkpoint) >0){
            print(volcanoplot(
            efit.15,
            coef      = coef,
            style     = "p-value",
            main      = coef,
            xlab      = "Log fold change",
            ylab      = "Log odds"
            ))
            lp  <- -log10(efit.15$p.value[row.names(results.15[which(results.15[,coef]!=0),]),coef])
            ord <- order(lp, decreasing = TRUE)
            x   <- efit.15$coef[c(names(lp[ord])),coef]
            y   <- lp[ord] + 0.1
            points(efit.15$coef[ord,coef], lp[ord], pch = 16, cex = 0.6, col = 2)
            text(x, y, 
                labels = paste(row.names(efit.15)[ord], " (",
                                efit.15$genes$COMMON[match(row.names(efit.15)[ord], efit.15$genes$ENSEMBL)], ")",
                                sep = ""), cex=0.4, col="red")
        } else{
            print(volcanoplot(
            efit.15,
            coef      = coef,
            style     = "p-value",
            main      = coef,
            xlab      = "Log fold change",
            ylab      = "Log odds"
            ))
        }
    }
dev.off()

# toptables for each coef + expression plots

for (coef in colnames(efit.15)) {
    tt <- topTable(
        efit.15,
        coef = coef,
        adjust.method = 'BH',
        sort.by = "M",
        resort.by = "B",
        p.value = 0.1,
        lfc = log2(1),
        number = Inf
    )
    if (length(tt) > 0)    {
        write.csv(
            tt, 
            file = paste("./toptable_", coef, "_15_p1.csv", sep =""),
            row.names = TRUE
        )
        plot_expression_tt(v.15, tt, "expgroup", "./15_p1_", coef, page_height = 20)

    }
}

### now plotting the DE results from the decideTests object
coef <- "E_L_15vsC_15"
pdf("./heatmap_E_L_15_vs-C.pdf", height = 7)
    prettyheat_15(c(row.names(results.15[which(results.15[,coef]<0),]), row.names(results.15[which(results.15[,coef]>0),])),
                            v.15,
                            metadata,
                            efit.15,
                            "E-L vs Control")
#dev.off()

coef <- "E_R_15vsC_15"
pdf("./heatmap_E_R_15_vs-C.pdf", height = 3)
    prettyheat_15(c(row.names(results.15[which(results.15[,coef]<0),]), row.names(results.15[which(results.15[,coef]>0),])),
                            v.15,
                            metadata,
                            efit.15,
                            "E-R vs Control")
#dev.off()


############################################################
# --------------------- Stage: E 17.5 -----------------------
############################################################

dir.create("/fs03/ag36/Shani/Trev-Seq/featurecounts_edgeR_Mar_2023/E17", showWarnings = FALSE)

setwd("/fs03/ag36/Shani/Trev-Seq/featurecounts_edgeR_Mar_2023/E17")

dge.17 <- dge[,as.character(metadata$ID.NEW[metadata$STAGE == "E17_5"])]

# filtering
dge.17.full  <- dge.17 # keep the old one in case we mess up
head(dge.17$counts)
head(cpm(dge.17))

# filter
keepMin <- apply(dge.17$counts, 1, max) >= 0.0
keepCpm <- rowSums(cpm(dge.17)>1) >= 2                  # Keep only genes with cpm above x in at least y samples
keep    <- keepMin & keepCpm
dge.17  <- dge.17[keep,]

dim(dge.17)
# [1] 20102    10

dge.17$samples$lib.size <- colSums(dge.17$counts)

# TMM normalisation
dge.17   <- calcNormFactors(dge.17)
dge.17

# normalisation factors
tmp.df           <- as.data.frame(cbind(dge.17$samples$norm.factors,rownames(dge.17$samples))) 
colnames(tmp.df) <- c("norm.factors", "sample")
tmp.df$month     <- as.factor(substr(tmp.df$sample, 1,3))
lcpm             <- cpm(dge.17, log=TRUE)

pdf("./norm.factors.normalised.17.CC.pdf")
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

pdf("./norm.factors.boxplot.normalised.17.CC.pdf")
   colourCount  <- length(unique(dge.17$samples$genotype.s))
   getPalette   <- colorRampPalette(brewer.pal(7, "Dark2"))
   col_groups   <- getPalette(colourCount)

   par(cex.axis=0.45)
   boxplot(lcpm, las=2, col=col_groups[dge.17$samples$genotype.s], main="")
   title(ylab="Log-cpm")
dev.off()

# MDS plot before voom fitting
pdf("./mds.plot.17.CC.pdf", width = 10, height = 10)
    plotMDS(dge.17, 
            labels = dge.17$samples$m.samples.s,
            col=as.numeric(dge.17$samples$genotype),
            cex = 0.5)
    legend("bottomleft",
        legend = as.character(unique(dge.17$samples$genotype)),
        col = as.numeric(unique(dge.17$samples$genotype)),
        pch = 15,
        cex = 0.7,
        pt.cex = 0.5)
dev.off()

# design matrix
group.17            <- as.factor(metadata$expgroup[match(colnames(dge.17$counts), metadata$ID.NEW)])
design.17           <- model.matrix(~0+group.17)
colnames(design.17) <- gsub("group\\.17","", colnames(design.17))

# adding surrogate variables
sv <- wsva(cpm(dge.17, log=TRUE, prior.count=2),
        design=design.17, n.sv = 4, 
        weight.by.sd = TRUE, plot = FALSE)

design.17           <- model.matrix(~0+group.17 + sv)
colnames(design.17) <- gsub("group\\.17","", colnames(design.17))

# BCV plot
dge.17 <- estimateDisp(dge.17, design.17, robust=TRUE)
pdf("bcv.17.CC.pdf")
    plotBCV(dge.17)
dev.off()

# MDS plot before voom fitting
pdf("./mds.plot.17.final.CC.pdf", width = 10, height = 10)
    plotMDS(dge.17, 
            labels = dge.17$samples$m.samples.s,
            col=as.numeric(dge.17$samples$genotype),
            cex = 1.5)
    legend("bottomleft",
        legend = as.character(unique(dge.17$samples$genotype)),
        col = as.numeric(unique(dge.17$samples$genotype)),
        pch = 15,
        cex = 1.0,
        pt.cex = 1.5)
dev.off()

# voom
pdf("./voom.dge.17.normalised.CC.pdf")
  v.17   <- voomWithQualityWeights(dge.17, design.17, plot = TRUE)
dev.off()

# Fitting linear model
fit.17 <- lmFit(object=v.17, design=design.17)

#  contrast matrix 
cont.matrix.17 <- makeContrasts(
  # at E17.5 
  E_L_17vsE_R_17 = E_L_17-E_R_17, 
  E_L_17vsC_17   = E_L_17-C_17, 
  E_R_17vsC_17   = E_R_17-C_17,
  E_L_17vsRest   = E_L_17-(E_R_17+C_17),
  E_R_17vsRest   = E_R_17-(E_L_17+C_17),
  E_17VsC_17     = (E_L_17+E_R_17)-C_17,
  levels=colnames(design.17)) 

# contrast fitting
fit2.17 <- contrasts.fit(fit.17, cont.matrix.17)
efit.17 <- eBayes(fit2.17)

# look at DE results
results.17 <- decideTests(efit.17, 
                       method = "separate", 
                       adjust.method = "fdr",
                       p.value=0.05,
                       lfc=0)

print(summary(results.17))

#        E_L_17vsE_R_17 E_L_17vsC_17 E_R_17vsC_17 E_L_17vsRest E_R_17vsRest
# Down                0            0            3        12449        12663
# NotSig          20102        20102        20091         7624         7424
# Up                  0            0            8           29           15
#        E_17VsC_17
# Down           20
# NotSig       7497
# Up          12585

write.table(summary(results.17),
            file      = "./de_summary_17_p1.CC.txt",
            quote     = FALSE,
            sep       = "\t",
            row.names = TRUE)

pdf("./efit.17.CC.pdf")
  plotSA(efit.17)
dev.off()

# save the v.17, dge.17 and efit.17
saveRDS(dge.17, "./dge17.final.rds")
saveRDS(v.17, "./voom17.final.rds")
saveRDS(efit.17, "./efit17.final.rds")
saveRDS(design.17, "./design17.final.rds")

# top table only for the gene lists, etc
out.17            <- topTreat(efit.17, adjust="BH", n = Inf, sort.by = "logFC", conf=0.99)
write.table(out.17, "./toptable_17_CC_new_2.csv", sep=",", quote = FALSE, row.names = TRUE)

# MD plot
pdf("./MD_plots_all_p1_Nov_17.CC_toptreat.pdf")   
    for (coef in colnames(efit.17)) {
        checkpoint <- row.names(results.17[which(results.17[,coef]!=0),])
            if (length(checkpoint) >0){
            status <- ifelse(efit.17$genes$ENSEMBL %in% row.names(results.17[which(results.17[,coef]<0),]), "downregulated", 
                                ifelse(efit.17$genes$ENSEMBL %in% row.names(results.17[which(results.17[,coef]>0),]), "upregulated", "not-sig"))
            values <- c("downregulated", "upregulated", "not-sig")
            hl.col <- c("blue","red","grey90")
            plotMD(efit.17, 
                    coef = coef, 
                    main = coef, 
                    cex = c(0.7,0.7,0.3), 
                    status =status, values=values, hl.col=hl.col,
                    args.legend = list(x = "topright", bty = "n", inset=c(-0.15, 0)),
                    tag.genes = efit.17$genes$ENSEMBL %in% row.names(results.17[which(results.17[,coef]!=0),])
            )
        }
    }
dev.off()

# interactive plots - MD
for (coef.id in 1:length(colnames(efit.17))) {
    md_plot_interactive(efit.17, dge.17, coef.id, paste0(colnames(efit.17)[coef.id],"_Nov_CC_new")) 
}

# volcano plot
pdf("./volcano_plots_all_p1_17_Nov.CC_new.pdf")
    for (coef in colnames(efit.17)) {
        checkpoint <- row.names(results.17[which(results.17[,coef]!=0),])
        if (length(checkpoint) >0){
            print(volcanoplot(
            efit.17,
            coef      = coef,
            style     = "p-value",
            main      = coef,
            xlab      = "Log fold change",
            ylab      = "Log odds"
            ))
            lp  <- -log10(efit.17$p.value[row.names(results.17[which(results.17[,coef]!=0),]),coef])
            ord <- order(lp, decreasing = TRUE)
            x   <- efit.17$coef[c(names(lp[ord])),coef]
            y   <- lp[ord] + 0.1
            points(efit.17$coef[ord,coef], lp[ord], pch = 16, cex = 0.6, col = 2)
            text(x, y, 
                labels = paste(row.names(efit.17)[ord], " (",
                                efit.17$genes$COMMON[match(row.names(efit.17)[ord], efit.17$genes$ENSEMBL)], ")",
                                sep = ""), cex=0.4, col="red")
        } else{
            print(volcanoplot(
            efit.17,
            coef      = coef,
            style     = "p-value",
            main      = coef,
            xlab      = "Log fold change",
            ylab      = "Log odds"
            ))
        }
    }
dev.off()

# toptables for each coef + expression plots

for (coef in colnames(efit.17)) {
    tt <- topTable(
        efit.17,
        coef = coef,
        adjust.method = 'BH',
        sort.by = "M",
        resort.by = "B",
        p.value = 0.1,
        lfc = log2(1),
        number = Inf
    )
    if (length(tt) > 0)    {
        write.csv(
            tt, 
            file = paste("./toptable_", coef, "_17_p1.csv", sep =""),
            row.names = TRUE
        )
        plot_expression_tt(v.17, tt, "expgroup", "./17_p1_", coef, page_height = 20)

    }
}

# now plotting the DE results from the decideTests object
coef <- "E_L_17vsC_17"
pdf("./heatmap_E_L_17_vs-C.pdf", height = 5)
    prettyheat_17(c(row.names(results.17[which(results.17[,coef]<0),]), row.names(results.17[which(results.17[,coef]>0),])),
                            v.17,
                            metadata,
                            efit.17,
                            "E-L vs Control")
#dev.off()

coef <- "E_R_17vsC_17"
pdf("./heatmap_E_R_17_vs-C.pdf", height = 5)
    prettyheat_17(c(row.names(results.17[which(results.17[,coef]<0),]), row.names(results.17[which(results.17[,coef]>0),])),
                            v.17,
                            metadata,
                            efit.17,
                            "E-R vs Control")
#dev.off()

############################################################
# --------------------- Stage: P0 -----------------------
############################################################

dir.create("/fs03/ag36/Shani/Trev-Seq/featurecounts_edgeR_Mar_2023/P0", showWarnings = FALSE)
setwd("/fs03/ag36/Shani/Trev-Seq/featurecounts_edgeR_Mar_2023/P0")

dge.0 <- dge[,as.character(metadata$ID.NEW[metadata$STAGE == "P0"])]

# filtering
dge.0.full  <- dge.0 # keep the old one in case we mess up
head(dge.0$counts)
head(cpm(dge.0))

table(rowSums(dge.0$counts==0)==2)

# filter
group  <- as.factor(metadata$expgroup[match(colnames(dge.0$counts), metadata$ID.NEW)])

keep   <- filterByExpr(dge.0, group=group)
dge.0 <- dge.0[keep,]

dim(dge.0)
# [1] 19692    10

dge.0$samples$lib.size <- colSums(dge.0$counts)

# TMM normalisation
dge.0   <- calcNormFactors(dge.0)
dge.0

# normalisation plots
tmp.df           <- as.data.frame(cbind(dge.0$samples$norm.factors,rownames(dge.0$samples))) 
colnames(tmp.df) <- c("norm.factors", "sample")
tmp.df$month     <- as.factor(substr(tmp.df$sample, 1,3))
lcpm             <- cpm(dge.0, log=TRUE)

pdf("./norm.factors.normalised.0.CC.pdf")
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

pdf("./norm.factors.boxplot.normalised.0.CC.pdf")
   colourCount  <- length(unique(dge.0$samples$genotype.s))
   getPalette   <- colorRampPalette(brewer.pal(7, "Dark2"))
   col_groups   <- getPalette(colourCount)

   par(cex.axis=0.45)
   boxplot(lcpm, las=2, col=col_groups[dge.0$samples$genotype.s], main="")
   title(ylab="Log-cpm")
dev.off()

# MDS plot
pdf("./mds.plot.0.CC.pdf", width = 10, height = 10)
    plotMDS(dge.0, 
            labels = dge.0$samples$m.samples.s,
            col=as.numeric(dge.0$samples$genotype),
            cex = 0.5)
    legend("bottomleft",
        legend = as.character(unique(dge.0$samples$genotype)),
        col = as.numeric(unique(dge.0$samples$genotype)),
        pch = 15,
        cex = 0.7,
        pt.cex = 0.5)
dev.off()

# design matrix
group.0            <- as.factor(metadata$expgroup[match(colnames(dge.0$counts), metadata$ID.NEW)])
design.0           <- model.matrix(~0+group.0)
colnames(design.0) <- gsub("group\\.0","", colnames(design.0))

# adding surrogate variables
sv <- wsva(cpm(dge.0, log=TRUE, prior.count=2),
        design=design.0, n.sv = 1, 
        weight.by.sd = TRUE, plot = FALSE)

design.0           <- model.matrix(~0+group.0 + sv)
colnames(design.0) <- gsub("group\\.0","", colnames(design.0))

# MDS plot before voom fitting
pdf("./mds.plot.0.final.CC.pdf", width = 10, height = 10)
    plotMDS(dge.0, 
            labels = dge.0$samples$m.samples.s,
            col=as.numeric(dge.0$samples$genotype),
            cex = 0.5)
    legend("bottomleft",
        legend = as.character(unique(dge.0$samples$genotype)),
        col = as.numeric(unique(dge.0$samples$genotype)),
        pch = 15,
        cex = 0.7,
        pt.cex = 0.5)
dev.off()

# BCV plot
dge.0 <- estimateDisp(dge.0, design.0, robust=TRUE)
pdf("bcv.0.CC.pdf")
    plotBCV(dge.0)
dev.off()

# voom
pdf("./voom.dge.0.normalised.CC.pdf")
  v.0   <- voomWithQualityWeights(dge.0, design.0, plot = TRUE)
dev.off()

# Fitting linear model
fit.0 <- lmFit(object=v.0, design=design.0)

#  contrast matrix 
cont.matrix.0 <- makeContrasts(
  # at E0.5 
  E_L_0vsE_R_0 = E_L_0-E_R_0, 
  E_L_0vsC_0   = E_L_0-C_0, 
  E_R_0vsC_0   = E_R_0-C_0,
  E_L_0vsRest  = E_L_0-(E_R_0+C_0),
  E_R_0vsRest  = E_R_0-(E_L_0+C_0),
  E_0VsC_0     = (E_L_0+E_R_0)-C_0,
  levels=colnames(design.0)) 

# contrast fitting
fit2.0 <- contrasts.fit(fit.0, cont.matrix.0)
efit.0 <- eBayes(fit2.0)

# look at DE results
results.0 <- decideTests(efit.0, 
                       method = "separate", 
                       adjust.method = "fdr",
                       p.value=0.05,
                       lfc=0)

print(summary(results.0))

#        E_L_0vsE_R_0 E_L_0vsC_0 E_R_0vsC_0 E_L_0vsRest E_R_0vsRest E_0VsC_0
# Down              0          0          0       13679       13766       10
# NotSig        19692      19692      19692        5956        5900     5873
# Up                0          0          0          57          26    13809

write.table(summary(results.0),
            file      = "./de_summary_0_p1.CC.txt",
            quote     = FALSE,
            sep       = "\t",
            row.names = TRUE)

pdf("./efit.0.CC.pdf")
  plotSA(efit.0)
dev.off()

# save the v.0, dge.0 and efit.0
saveRDS(dge.0, "./dge0.final.rds")
saveRDS(v.0, "./voom0.final.rds")
saveRDS(efit.0, "./efit0.final.rds")
saveRDS(design.0, "./design0.final.rds")

# top tables
out.0           <- topTable(efit.0, adjust="BH", p.value=0.1, lfc=log2(1), n = Inf)
write.table(out.0, "./toptable_0_p1.CC.csv", sep=",", quote = FALSE, row.names = TRUE)

# MD plot
pdf("./MD_plots_all_p1_0.CC.pdf")
    for (coef in colnames(efit.0)) {
        plotMD(efit.0, 
                coef = coef, 
                main = coef, 
                cex = c(0.3,0.5,0.5), 
                tag.genes = efit.0$genes$ENSEMBL %in% row.names(results.0[which(results.0[,coef]!=0),])
        )
    }
dev.off()

# interactive plots - MD
for (coef.id in 1:length(colnames(efit.0))) {
    md_plot_interactive(efit.0, dge.0, coef.id, paste0(colnames(efit.0)[coef.id], "_CC")) 
}

# volcano plot
pdf("./volcano_plots_all_p1_0.CC.pdf")
    for (coef in colnames(efit.0)) {
        checkpoint <- row.names(results.0[which(results.0[,coef]!=0),])
        if (length(checkpoint) >0){
            print(volcanoplot(
            efit.0,
            coef      = coef,
            style     = "p-value",
            main      = coef,
            xlab      = "Log fold change",
            ylab      = "Log odds"
            ))
            lp  <- -log10(efit.0$p.value[row.names(results.0[which(results.0[,coef]!=0),]),coef])
            ord <- order(lp, decreasing = TRUE)
            x   <- efit.0$coef[c(names(lp[ord])),coef]
            y   <- lp[ord] + 0.1
            points(efit.0$coef[ord,coef], lp[ord], pch = 16, cex = 0.6, col = 2)
            text(x, y, 
                labels = paste(row.names(efit.0)[ord], " (",
                                efit.0$genes$COMMON[match(row.names(efit.0)[ord], efit.0$genes$ENSEMBL)], ")",
                                sep = ""), cex=0.4, col="red")
        } else{
            print(volcanoplot(
            efit.0,
            coef      = coef,
            style     = "p-value",
            main      = coef,
            xlab      = "Log fold change",
            ylab      = "Log odds"
            ))
        }
    }
dev.off()

# toptables for each coef + expression plots
for (coef in colnames(efit.0)) {
    tt <- topTable(
        efit.0,
        coef = coef,
        adjust.method = 'BH',
        sort.by = "M",
        resort.by = "B",
        p.value = 0.1,
        lfc = log2(1),
        number = Inf
    )
    if (length(tt) > 0)    {
        write.csv(
            tt, 
            file = paste("./toptable_", coef, "_0_p1.CC.csv", sep =""),
            row.names = TRUE
        )
        plot_expression_tt(v.0, tt, "expgroup", "./0_p1_CC_", coef, 50)

    }
}

### now plotting the DE results from the decideTests object
coef <- "E_L_0vsRest"
tt.L <- topTable(
        efit.0,
        coef = coef,
        adjust.method = 'BH',
        sort.by = "M",
        resort.by = "B",
        p.value = 0.1,
        lfc = log2(1),
        number = Inf
    )
pdf("./heatmap_E_L_0_vs-Rest.pdf", height = 5)
    #prettyheat_0(c(row.names(results.0[which(results.0[,coef]<0),]), row.names(results.0[which(results.0[,coef]>0),]))[1:20],
    prettyheat_0(c(row.names(tt.L))[1:20],
                            v.0,
                            metadata,
                            efit.0,
                            "E-L vs Rest | Top 20")
#dev.off()
coef <- "E_R_0vsRest"
tt.R <- topTable(
        efit.0,
        coef = coef,
        adjust.method = 'BH',
        sort.by = "M",
        resort.by = "B",
        p.value = 0.1,
        lfc = log2(1),
        number = Inf
    )
pdf("./heatmap_E_R_0_vs-Rest.pdf", height = 9)
    #prettyheat_0(c(row.names(results.0[which(results.0[,coef]<0),]), row.names(results.0[which(results.0[,coef]>0),]))[1:20],
    prettyheat_0(c(row.names(tt.R))[1:20],
                            v.0,
                            metadata,
                            efit.0,
                            "E-R vs Rest | Top 20")
#dev.off()

############################################################
# --------------------- Stage: P3 -----------------------
############################################################

dir.create("/fs03/ag36/Shani/Trev-Seq/featurecounts_edgeR_Mar_2023/P3", showWarnings = FALSE)
setwd("/fs03/ag36/Shani/Trev-Seq/featurecounts_edgeR_Mar_2023/P3")

dge.3 <- dge[,as.character(metadata$ID.NEW[metadata$STAGE == "P3"])]

dge.3.full  <- dge.3 # keep the old one in case we mess up
head(dge.3$counts)
head(cpm(dge.3))

#apply(dge.3$counts, 1, max) # total gene counts per sample

table(rowSums(dge.3$counts==0)==2)

# filter
group  <- as.factor(metadata$expgroup[match(colnames(dge.3$counts), metadata$ID.NEW)])

keep   <- filterByExpr(dge.3, group=group)
dge.3 <- dge.3[keep,]

dim(dge.3)
# 

dge.3$samples$lib.size <- colSums(dge.3$counts)

# assigning gene common name data frame
dge.3$genes <- as.data.frame(cbind(ENSEMBL = row.names(dge.3),
                                   COMMON = genes$SYMBOL[match(row.names(dge.3), genes$ENSEMBL)]))

# TMM normalisation
dge.3   <- calcNormFactors(dge.3)
dge.3

# normalisation plots
tmp.df           <- as.data.frame(cbind(dge.3$samples$norm.factors,rownames(dge.3$samples))) 
colnames(tmp.df) <- c("norm.factors", "sample")
tmp.df$month     <- as.factor(substr(tmp.df$sample, 1,3))
lcpm             <- cpm(dge.3, log=TRUE)

pdf("./norm.factors.normalised.3.CC.pdf")
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

pdf("./norm.factors.boxplot.normalised.3.CC.pdf")
   colourCount  <- length(unique(dge.3$samples$genotype.s))
   getPalette   <- colorRampPalette(brewer.pal(7, "Dark2"))
   col_groups   <- getPalette(colourCount)

   par(cex.axis=0.45)
   boxplot(lcpm, las=2, col=col_groups[dge.3$samples$genotype.s], main="")
   title(ylab="Log-cpm")
dev.off()

# MDS plot
pdf("./mds.plot.3.CC.pdf", width = 10, height = 10)
    plotMDS(dge.3, 
            labels = dge.3$samples$m.samples.s,
            col=as.numeric(dge.3$samples$genotype),
            cex = 0.5)
    legend("bottomleft",
        legend = as.character(unique(dge.3$samples$genotype)),
        col = as.numeric(unique(dge.3$samples$genotype)),
        pch = 15,
        cex = 0.7,
        pt.cex = 0.5)
dev.off()

# design matrix
group.3            <- as.factor(metadata$expgroup[match(colnames(dge.3$counts), metadata$ID.NEW)])
design.3           <- model.matrix(~0+group.3)
colnames(design.3) <- gsub("group\\.3","", colnames(design.3))

# adding surrogate variables
sv <- wsva(cpm(dge.3, log=TRUE, prior.count=2),
        design=design.3, n.sv = 1, 
        weight.by.sd = TRUE, plot = FALSE)

design.3           <- model.matrix(~0+group.3 + sv)
colnames(design.3) <- gsub("group\\.3","", colnames(design.3))

# BCV plot
dge.3 <- estimateDisp(dge.3, design.3, robust=TRUE)
pdf("bcv.3.CC.pdf")
    plotBCV(dge.3)
dev.off()

# MDS plot before voom fitting
pdf("./mds.plot.3.final.CC.pdf", width = 10, height = 10)
    plotMDS(dge.3, 
            labels = dge.3$samples$m.samples.s,
            col=as.numeric(dge.3$samples$genotype),
            cex = 0.5)
    legend("bottomleft",
        legend = as.character(unique(dge.3$samples$genotype)),
        col = as.numeric(unique(dge.3$samples$genotype)),
        pch = 15,
        cex = 0.7,
        pt.cex = 0.5)
dev.off()

# voom
pdf("./voom.dge.3.normalised.CC.pdf")
    v.3   <- voomWithQualityWeights(dge.3, design.3, plot = TRUE)
dev.off()

# Fitting linear model
fit.3 <- lmFit(object=v.3, design=design.3)

#  contrast matrix 
cont.matrix.3 <- makeContrasts(
  # at E0.5 
  E_L_3vsE_R_3 = E_L_3-E_R_3, 
  E_L_3vsC_3   = E_L_3-C_3, 
  E_R_3vsC_3   = E_R_3-C_3,
  E_L_3vsRest  = E_L_3-(E_R_3+C_3),
  E_R_3vsRest  = E_R_3-(E_L_3+C_3),
  E_3VsC_3     = (E_L_3+E_R_3)-C_3,
  levels=colnames(design.3)) 

# contrast fitting
fit2.3 <- contrasts.fit(fit.3, cont.matrix.3)
efit.3 <- eBayes(fit2.3)

# look at DE results
results.3 <- decideTests(efit.3, 
                       method = "separate", 
                       adjust.method = "fdr",
                       p.value=0.05,
                       lfc=0)

print(summary(results.3))

#           E_L_3vsE_R_3 E_L_3vsC_3 E_R_3vsC_3 E_L_3vsRest E_R_3vsRest E_3VsC_3
# Down              0          0          0       13364       13511       20
# NotSig        19651      19649      19644        6224        6078     6196
# Up                0          2          7          63          62    13435


write.table(summary(results.3),
            file      = "./de_summary_3_p1.CC.txt",
            quote     = FALSE,
            sep       = "\t",
            row.names = TRUE)

pdf("./efit.3.CC.pdf")
  plotSA(efit.3)
dev.off()

# save the v.3, dge.3 and efit.3
saveRDS(dge.3, "./dge3.final.rds")
saveRDS(v.3, "./voom3.final.rds")
saveRDS(efit.3, "./efit3.final.rds")
saveRDS(design.3, "./design3.final.rds")

# top tables
out.3           <- topTable(efit.3, adjust="BH", p.value=0.1, lfc=log2(1), n = Inf)
write.table(out.3, "./toptable_3_p1.CC.csv", sep=",", quote = FALSE, row.names = TRUE)

# MD plot
pdf("./MD_plots_all_p1_3.CC.pdf")
    for (coef in colnames(efit.3)) {
        plotMD(efit.3, 
                coef = coef, 
                main = coef, 
                cex = c(0.3,0.5,0.5), 
                tag.genes = efit.3$genes$ENSEMBL %in% row.names(results.3[which(results.3[,coef]!=0),])
        )
    }
dev.off()

# interactive plots - MD
for (coef.id in 1:length(colnames(efit.3))) {
    md_plot_interactive(efit.3, dge.3, coef.id, paste0(colnames(efit.3)[coef.id], "_CC")) 
}

# volcano plot
pdf("./volcano_plots_all_p1_3.CC.pdf")
    for (coef in colnames(efit.3)) {
        checkpoint <- row.names(results.3[which(results.3[,coef]!=0),])
        if (length(checkpoint) >0){
            print(volcanoplot(
            efit.3,
            coef      = coef,
            style     = "p-value",
            main      = coef,
            xlab      = "Log fold change",
            ylab      = "Log odds"
            ))
            lp  <- -log10(efit.3$p.value[row.names(results.3[which(results.3[,coef]!=0),]),coef])
            ord <- order(lp, decreasing = TRUE)
            x   <- efit.3$coef[c(names(lp[ord])),coef]
            y   <- lp[ord] + 0.1
            points(efit.3$coef[ord,coef], lp[ord], pch = 16, cex = 0.6, col = 2)
            text(x, y, 
                labels = paste(row.names(efit.3)[ord], " (",
                                efit.3$genes$COMMON[match(row.names(efit.3)[ord], efit.3$genes$ENSEMBL)], ")",
                                sep = ""), cex=0.4, col="red")
        } else{
            print(volcanoplot(
            efit.3,
            coef      = coef,
            style     = "p-value",
            main      = coef,
            xlab      = "Log fold change",
            ylab      = "Log odds"
            ))
        }
    }
dev.off()

# toptables for each coef + expression plots

for (coef in colnames(efit.3)) {
    tt <- topTable(
        efit.3,
        coef = coef,
        adjust.method = 'BH',
        sort.by = "M",
        resort.by = "B",
        p.value = 0.1,
        lfc = log2(1),
        number = Inf
    )
    if (length(tt) > 0)    {
        write.csv(
            tt, 
            file = paste("./toptable_", coef, "_3_p1.CC.csv", sep =""),
            row.names = TRUE
        )
        plot_expression_tt(v.3, tt, "expgroup", "./3_p1_CC_", coef, 50)

    }
}

### now plotting the DE results from the decideTests object
coef <- "E_L_3vsRest"
tt.L <- topTable(
        efit.3,
        coef = coef,
        adjust.method = 'BH',
        sort.by = "M",
        resort.by = "B",
        p.value = 0.1,
        lfc = log2(1),
        number = Inf
    )
pdf("./heatmap_E_L_3_vs-Rest.pdf", height = 5)
    #prettyheat_3(c(row.names(results.3[which(results.3[,coef]<0),]), row.names(results.3[which(results.3[,coef]>0),]))[1:20],
    prettyheat_3(c(row.names(tt.L))[1:20],
                            v.3,
                            metadata,
                            efit.3,
                            "E-L vs Rest | Top 20")
#dev.off()
coef <- "E_R_3vsRest"
tt.R <- topTable(
        efit.3,
        coef = coef,
        adjust.method = 'BH',
        sort.by = "M",
        resort.by = "B",
        p.value = 0.1,
        lfc = log2(1),
        number = Inf
    )
pdf("./heatmap_E_R_3_vs-Rest.pdf", height = 9)
    #prettyheat_3(c(row.names(results.3[which(results.3[,coef]<0),]), row.names(results.3[which(results.3[,coef]>0),]))[1:20],
    prettyheat_3(c(row.names(tt.R))[1:20],
                            v.3,
                            metadata,
                            efit.3,
                            "E-R vs Rest | Top 20")
#dev.off()
