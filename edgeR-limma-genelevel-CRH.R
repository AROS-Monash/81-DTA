#' ----- Ecrh.dgeR-Limma pipeline on Trev-Seq dataset (gene level) --------
#' @author Shani Amaraisnghe
#' @date 08/08/2022
#' @aliases Rosello-Diez lab, ARMI, Monash
#' @concept analysing Trev-seq bulk RNA-Seq data for Chee
#' @description testing the similarities between the Adelaide CRH dataset
#' and the Trev-Seq dataset
#' Received the count matrix and metadata from Vincent
#' running the ecrh.dgeR-limma pipeline in M3
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
library(Glimma)
library(RColorBrewer)
library(Mus.musculus)
library(plyr) 

# -------------I/O --------------
# working in M3
#module load R/4.2-mkl

# ------------- Set working directory ----

# location of processing
setwd("/fs03/ag36/Shani/Trev-Seq/CRH")
here::here()

# ------------ Read in raw count file ------------

crh.counts  <- read.table("/fs03/ag36/Shani/Trev-Seq/CRH/rawCounts.txt", header = TRUE)

# ------------- Read in the metadata file -----------------

crh.metadata            <- read.table("/fs03/ag36/Shani/Trev-Seq/CRH/sampleInfo.csv", sep = ",", header = TRUE)
crh.metadata$sampleName <- str_replace(crh.metadata$sampleName, "-", "\\.")
crh.metadata$label      <- factor(crh.metadata$label)

# -------------- Adjust the column names for creating the crh.dgelist object -----------

crh.dge  <- DGEList(counts = crh.counts[, -1], 
                    group  = as.factor(crh.metadata$label[match(colnames(crh.counts[, -1]), crh.metadata$sampleName)]),
                    gene   = crh.counts[,1])

rownames(crh.dge)         <- crh.dge$genes$genes
crh.dge$samples$label     <- crh.metadata$label[match(colnames(crh.dge), crh.metadata$sampleName)]
crh.dge$samples$group     <- crh.metadata$group[match(colnames(crh.dge), crh.metadata$sampleName)]
crh.dge$samples$treatment <- crh.metadata$treatment[match(colnames(crh.dge), crh.metadata$sampleName)]
crh.dge$samples$mouseID   <- crh.metadata$mouseID[match(colnames(crh.dge), crh.metadata$sampleName)]
crh.dge$samples$timepoint <- crh.metadata$timepoint[match(colnames(crh.dge), crh.metadata$sampleName)]

# --------- Plot library size before normalisation ----------

tmp.df           <- as.data.frame(cbind(lib.size     = crh.dge$samples$lib.size,
                                        labels       = crh.dge$samples$labels, 
                                        group        = crh.dge$samples$group,
                                        mouseID      = crh.dge$samples$mouseID, 
                                        timepoint    = crh.dge$samples$timepoint, 
                                        treatment    = crh.dge$samples$treatment))

tmp.df$id.treatment <- paste(tmp.df$mouseID, tmp.df$treatment, sep = "-")                                         
  
pdf("./library.size.unnormalised.full.crh.pdf")
  ggplot(tmp.df, aes(x    = as.factor(id.treatment), 
                     y    = as.numeric(lib.size),
                     fill = as.factor(group))) + 
    geom_bar(stat = "identity") +
    theme(legend.position="none") +
    scale_fill_manual(values = c("deepskyblue3","darkorange3","goldenrod1", "deepskyblue4","darkorange4","goldenrod4")) +
    facet_grid(~as.factor(timepoint), scales="free", space="free_x",) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ylab("library size")
dev.off()

pdf("./log.library.size.unnormalised.full.crh.pdf")
ggplot(tmp.df, aes(x    = as.factor(id.treatment), 
                   y    = log(as.numeric(lib.size)),
                   fill = as.factor(group))) + 
  geom_bar(stat = "identity") +
  theme(legend.position="none") +
  scale_fill_manual(values = c("deepskyblue3","darkorange3","goldenrod1", "deepskyblue4","darkorange4","goldenrod4")) +
  facet_grid(~as.factor(timepoint), scales="free", space="free_x",) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("log(library size)")
dev.off()

# ----------- MDS plot before normalisation -------------

pdf("./mds.plot.unnormalised.full.crh.pdf", width = 10, height = 10)
plotMDS(crh.dge,
        labels = crh.dge$samples$label,
        method ="logFC", 
        col    = as.numeric(as.factor(crh.dge$samples$timepoint)), 
        cex    = 0.5)
legend("bottomleft",
       legend = as.character(unique(crh.dge$samples$timepoint)),
       col    = as.numeric(unique(as.factor(crh.dge$samples$timepoint))),
       pch    = 15,
       cex    = 0.7,
       pt.cex = 0.5)
dev.off()

# colouring based on treatment and timepoint
pdf("./mds.plot.unnormalised.tt.full.crh.pdf", width = 10, height = 10)
plotMDS(crh.dge,
        labels = as.character(crh.dge$samples$label),
        method ="logFC", 
        col    = as.numeric(as.factor(paste(crh.dge$samples$treatment,crh.dge$samples$timepoint, sep = "-"))), 
        cex    = 0.5)
legend("bottomleft",
       legend = as.character(unique(as.factor(paste(crh.dge$samples$treatment,crh.dge$samples$timepoint, sep = "-")))),
       col    = as.numeric(unique(as.factor(paste(crh.dge$samples$treatment,crh.dge$samples$timepoint, sep = "-")))),
       pch    = 15,
       cex    = 0.7,
       pt.cex = 0.5)
dev.off()

# getting a sense of the sizes in the crh.dge object
L <- mean(crh.dge$samples$lib.size) * 1e-6
M <- median(crh.dge$samples$lib.size) * 1e-6
c(L, M)
# [1] 19.58112 19.79210


# --------------- Filter the crh.dgelist object based on transcript counts ----------------
dim(crh.dge)
# [1] 28467    36
crh.dge.full             <- crh.dge # keep the old one in case we mess up
head(crh.dge$counts)
head(cpm(crh.dge))

#apply(crh.dge$counts, 2, sum) # total gene counts per sample

table(rowSums(crh.dge$counts==0)==2)

# FALSE  TRUE 
# 28139   328 

# filter option 1
# keepMin <- apply(crh.dge$counts, 1, max) >= 0.0
# keepCpm <- rowSums(cpm(crh.dge)> 0.0) >= 2                  # Keep only genes with cpm above x in at least y samples
# keep    <- keepMin & keepCpm
# crh.dge     <- crh.dge[keep,]

# filter option 2 - which filters out all the DE genes on the final dataset

group  <- as.factor(crh.metadata$label[match(colnames(crh.dge$counts), crh.metadata$sampleName)])

keep <- filterByExpr(crh.dge, group=group)
crh.dge  <- crh.dge[keep,]

dim(crh.dge)
# [1] 16265    36

crh.dge$samples$lib.size <- colSums(crh.dge$counts)

# getting a sense of the sizes in the edited crh.dge object
L <- mean(crh.dge$samples$lib.size) * 1e-6
M <- median(crh.dge$samples$lib.size) * 1e-6
c(L, M)
# [1] 19.57440 19.78514

# ----------------- Assesing the batch effect for the total dataset -----------------

# Create design matrix for batch effect detection

group      <- as.factor(crh.metadata$label[match(colnames(crh.dge$counts), crh.metadata$sampleName)])
timepoint  <- as.factor(crh.metadata$timepoint[match(colnames(crh.dge$counts), crh.metadata$timepoint)])

crh.design <- model.matrix(~0+group) 

# run removebatcheffects()

tmp.crh.dge  <- crh.dge
tmp.crh.dge  <- calcNormFactors(tmp.crh.dge)
tmp.crh.dge

crh.dge.logcpm <- cpm(tmp.crh.dge,log=TRUE,prior.count=2)
crh.dge.logcpm <- removeBatchEffect(crh.dge.logcpm, batch=as.character(tmp.crh.dge$samples$timepoint), design = crh.design)

# MDS plot by timepoint
pdf("./mds.plot.batcheffectcorr.crh.pdf", width = 10, height = 10)
    plotMDS(crh.dge.logcpm, 
            labels = crh.dge$samples$label,
            col    = as.numeric(as.factor(crh.dge$samples$timepoint)), 
            cex    = 0.5)
    legend("bottomleft",
        legend = as.character(unique(crh.dge$samples$timepoint)),
        col    = as.numeric(unique(as.factor(crh.dge$samples$timepoint))),
        pch    = 15,
        cex    = 0.7,
        pt.cex = 0.5)
dev.off()

# Does not change anything - so no batch effect is corrected here

# -------------- TMM normalisation -----------------------

# The calcNormFactors() function normalizes for RNA composition by finding a set 
# of scaling factors for the library sizes that minimize the log-fold changes 
# between the samples for most genes. The default method for computing these 
# scale factors uses a trimmed mean of M-values (TMM) between each pair of 
# samples. We call the product of the original library size and the scaling 
# factor the effective library size.The effective library size replaces 
# the original library size in all downsteam analyses.
crh.dge  <- calcNormFactors(crh.dge)
crh.dge

# -------------- Extracting common gene names from Ensembl for annotations of DE genes --------

genes <- select(Mus.musculus, keys=row.names(crh.dge), columns=c("SYMBOL", "TXCHROM", "ENSEMBL"), 
                keytype="ENSEMBL")
dim(genes)

crh.dge$genes$common_name <- genes$SYMBOL[match(row.names(crh.dge), genes$ENSEMBL)]

saveRDS(crh.dge, "./crh.dge.full.final.rds")

# -------------- MDS plot for the normalised samples ------------

pdf("./mds.plot.normalised.crh.pdf", width = 10, height = 10)
plotMDS(crh.dge,
        labels = crh.dge$samples$samples.s,
        method ="logFC", 
        col    = as.numeric(as.factor(crh.dge$samples$timepoint)), 
        cex    = 0.5)
    legend("bottomleft",
        legend = as.character(unique(crh.dge$samples$timepoint)),
        col    = as.numeric(unique(as.factor(crh.dge$samples$timepoint))),
        pch    = 15,
        cex    = 0.7,
        pt.cex = 0.5)
dev.off()

# colouring based on treatment and timepoint

pdf("./mds.plot.tt.normalised.crh.pdf", width = 10, height = 10)
   colourCount  <- length(unique(as.numeric(unique(as.factor(paste(crh.dge$samples$treatment,crh.dge$samples$timepoint, sep = "-"))))))
   getPalette   <- colorRampPalette(brewer.pal(7, "Dark2"))
   col_groups   <- getPalette(colourCount)

plotMDS(crh.dge,
        labels = crh.dge$samples$label,
        method ="logFC", 
        col    = col_groups[as.numeric(as.factor(paste(crh.dge$samples$treatment,crh.dge$samples$timepoint, sep = "-")))],
        cex    = 0.5)
    legend("bottomleft",
        legend = as.character(unique(as.factor(paste(crh.dge$samples$treatment,crh.dge$samples$timepoint, sep = "-")))),
        col    = col_groups[as.numeric(as.factor(paste(crh.dge$samples$treatment,crh.dge$samples$timepoint, sep = "-")))],
        pch    = 15,
        cex    = 0.7,
        pt.cex = 0.5)
dev.off()

# ----------- Plotting the normalisation factors for final dataset ----------

tmp.df           <- as.data.frame(cbind(norm.factors = crh.dge$samples$norm.factors,
                                        labels       = crh.dge$samples$labels, 
                                        group        = crh.dge$samples$group,
                                        mouseID      = crh.dge$samples$mouseID, 
                                        timepoint    = crh.dge$samples$timepoint, 
                                        treatment    = crh.dge$samples$treatment))

tmp.df$id.treatment <- paste(tmp.df$mouseID, tmp.df$treatment, sep = "-")                                         
  
pdf("./norm.factors.normalised.full.crh.pdf")
    ggplot(tmp.df, aes(x    = as.factor(id.treatment), 
                    y    = as.numeric(norm.factors),
                    fill = as.factor(group))) + 
    geom_bar(stat = "identity") +
    theme(legend.position="none") +
    scale_fill_manual(values = c("deepskyblue3","darkorange3","goldenrod1", "deepskyblue4","darkorange4","goldenrod4")) +
    facet_grid(~as.factor(timepoint), scales="free", space="free_x",) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ylab("norm.factors")
dev.off()

lcpm <- cpm(crh.dge,log=TRUE,prior.count=2)

pdf("./norm.factors.boxplot.normalised.crh.pdf")
   colourCount  <- length(unique(as.numeric(unique(as.factor(paste(crh.dge$samples$treatment,crh.dge$samples$timepoint, sep = "-"))))))
   getPalette   <- colorRampPalette(brewer.pal(7, "Dark2"))
   col_groups   <- getPalette(colourCount)

   par(cex.axis=0.45)
   boxplot(lcpm, 
           las  = 2, 
           col  = col_groups[as.numeric(as.factor(paste(crh.dge$samples$treatment,crh.dge$samples$timepoint, sep = "-")))], 
           main = "")
   title(ylab = "Log-cpm")
dev.off()

# ----------------- Density plot for the counts --------------------

# here tpm.crh.dge is a normalised (TMM) crh.dge object of 56 samples (excludiing the sheath and the Dahups)
lcpm <- cpm(crh.dge,log=TRUE,prior.count=2)

pdf("./density.filtered.normalised.logcpm.crh.pdf")

    lcpm.cutoff <- log2(10/M + 2/L)
    nsamples    <- ncol(crh.dge)

    colourCount  <- length(unique(as.numeric(unique(as.factor(paste(crh.dge$samples$treatment,crh.dge$samples$timepoint, sep = "-"))))))    
    getPalette   <- colorRampPalette(brewer.pal(7, "Dark2"))
    col_groups   <- getPalette(colourCount)

    plot(density(lcpm[,1]), col=col_groups[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
    title(xlab="Log-cpm")
    #abline(v=lcpm.cutoff, lty=3)
    for (i in 2:nsamples){
      den <- density(lcpm[,i])
      lines(den$x, den$y, col=col_groups[i], lwd=2)
    }
    legend("topright", 
    legend = as.character(unique(colnames(tmp.crh.dge))),
    text.col=col_groups, bty="n",
    cex = 0.7)
dev.off()

# ----------------- save DG object and the metadata for future use ---------------

saveRDS(crh.dge, "/fs03/ag36/Shani/Trev-Seq/CRH/final.crh.dge.rds")
saveRDS(crh.metadata, "/fs03/ag36/Shani/Trev-Seq/CRH/final.metadata.rds")

# ----------------- Create design matrix for DE analysis ----------------------

# treatment  <- as.factor(crh.metadata$treatment[match(colnames(crh.dge$counts), crh.metadata$sampleName)])
# timepoint  <- as.factor(crh.metadata$timepoint[match(colnames(crh.dge$counts), crh.metadata$sampleName)])
# crh.design <- model.matrix(~0+treatment+timepoint)

mouse      <- as.factor(crh.metadata$sampleName[match(colnames(crh.dge$counts), crh.metadata$sampleName)])
group      <- as.factor(crh.metadata$group[match(colnames(crh.dge$counts), crh.metadata$sampleName)])
crh.design <- model.matrix(~0+group)

# ----------------- Create BCV plot ------------------------------------------

crh.dge <- estimateDisp(crh.dge, crh.design, robust=TRUE)
pdf("bcv.crh.pdf")
    plotBCV(crh.dge)
dev.off()

# -------------------- voom ---------------------

pdf("./voom.crh.dge.normalised.pdf")
  v   <- voom(crh.dge, crh.design, plot = TRUE)
dev.off()

# --------------------- Fitting linear model --------------

cor <- duplicateCorrelation(v,crh.design, block = mouse)
# fit <- lmFit(object=v, design=design, 
#  correlation=cor$consensus.correlation)

fit <- lmFit(object=v, design=crh.design)

# ------------------ contrast matrix --------------

cont.matrix <- makeContrasts(
  # at 24 hours 
  Unx.24vsSham.24 = groupUnx.24hr-groupSham.24hr, 

  # at 48 hours
  Unx.48vsSham.48 = groupUnx.48hr-groupSham.48hr, 
 
  # at 72 hours
  Unx.72vsSham.72 = groupUnx.72hr-groupSham.72hr, 

  # Unx vs Sham
  UnxvsSham       = (groupUnx.24hr + groupUnx.48hr + groupUnx.72hr) - (groupSham.24hr + groupSham.48hr + groupSham.72hr),

  # 24 hours vs 48 hours
  hr24vshr48      = (groupUnx.24hr-groupSham.24hr) - (groupUnx.48hr-groupSham.48hr),

  # 48 hours vs 72 hours
  hr48vshr72      = (groupUnx.48hr-groupSham.72hr) - (groupUnx.72hr-groupSham.72hr),

  # 24 hours vs 72 hours
  hr24vshr72      = (groupUnx.24hr-groupSham.24hr) - (groupUnx.72hr-groupSham.72hr),

  levels=colnames(crh.design)) 

# -------------------------- Contrast fitting --------------------

fit2 <- contrasts.fit(fit, cont.matrix)
efit <- eBayes(fit2)

results <- decideTests(efit, 
                       method = "separate", 
                       adjust.method = "fdr",
                       p.value=0.01,
                       lfc=1)

print(summary(results))

#        Unx.24vsSham.24 Unx.48vsSham.48 Unx.72vsSham.72 UnxvsSham hr24vshr48
# Down                40              95               3       224        235
# NotSig           16188           16034           16215     15582      15843
# Up                  37             136              47       459        187
#        hr48vshr72 hr24vshr72
# Down           48          6
# NotSig      16144      16257
# Up             73          2

write.table(summary(results),
            file      = "./de_summary_all_crh.txt",
            quote     = FALSE,
            sep       = "\t",
            row.names = TRUE)

pdf("./efit.filtered.crh.pdf")
  plotSA(efit)
dev.off()

# adding common gene names from Ensembl
genes <- select(Mus.musculus, keys=row.names(efit), columns=c("SYMBOL", "TXCHROM", "ENSEMBL"), 
                keytype="ENSEMBL")
dim(genes)

genes          <- genes[!is.na(genes$SYMBOL),]
efit$gene_name <- genes$SYMBOL[match(row.names(efit), genes$ENSEMBL)]

#out <- topTable(efit, n=Inf, adjust="BH")
out          <- topTable(efit, adjust="BH", p.value=0.01, lfc=1, n = Inf, sort.by = "F")
out$gene_name <- genes$SYMBOL[match(row.names(out), genes$ENSEMBL)]
out2 <- cbind(efit$coef[match(row.names(out), row.names(efit$coef))],
              out[, c('P.Value','adj.P.Val','AveExpr','F', 'gene_name')])

write.table(out2, "./toptable_all_crh.csv", sep=",", quote = FALSE, row.names = TRUE)

# ------------------------ Plot venn diagrams ---------------------------

pdf("./venn_interactions.crh.pdf", height = 11.69, width = 8.27, onefile = TRUE)
    # E15.5 contrast set
    print(vennDiagram(
        results[, c('Unx.24vsSham.24', 'Unx.48vsSham.48', 'Unx.72vsSham.72')],
        include = c("both"),
        counts.col = "black",
        cex = c(1.5,1,0.7) * 0.6,
        circle.col = c('skyblue', 'pink1', 'mediumorchid'),
        mar = rep(3,4),
        main = "Diff between timepoints"
    ))

    # E17.5 contrast set
    print(vennDiagram(
        results[, c('hr24vshr48', 'hr48vshr72', 'hr24vshr72')],
        include = c("both"),
        counts.col = "black",
        cex = c(1.5,1,0.7) * 0.6,
        circle.col = c('skyblue', 'pink1', 'mediumorchid'),
        mar = rep(3,4),
        main = "Diff between contrasts in timepoint 'interactions'"
    ))  
dev.off()

# ----------------- Plot volcano plots ----------------------

pdf("./volcano_plots_all_crh.pdf")
    for (coef in colnames(efit)) {
        print(volcanoplot(
            efit,
            coef      = coef,
            style     = "p-value",
            main      = coef,
            xlab      = "Log fold change",
            ylab      = "Log odds"
        ))
        lp  <- -log10(efit$p.value[,coef])
        ord <- order(lp, decreasing = TRUE)[1:10]
        x   <- efit$coef[ord,coef]
        y   <- lp[ord] + 0.7
        points(efit$coef[ord,coef], lp[ord], pch = 16, cex = 0.6, col = 2)
        text(x, y, efit$gene_name[ord], cex=0.4, col="red")
    }
dev.off()

# ----------------- Interactive plot - Volcano -----------------------------

for (coef.id in 1:length(colnames(efit))) {
    volcano_plot_interactive(efit[,coef.id], crh.dge, main = c(colnames(efit)[coef.id]), html = paste0("volcano_", c(colnames(efit)[coef.id]))) 
} #efit, dge, main = "", html = "volcano"

# ----------------- Interactive plot - MD -----------------------------

for (coef.id in 1:length(colnames(efit))) {
    main = c(colnames(efit)[coef.id]), html = paste0("volcano_", c(colnames(efit)[coef.id])
}

# ------------------- toptables for each coef + expression plots -------------------

for (coef in colnames(efit)) {
        tt <- topTable(
            efit,
            coef          = coef,
            adjust.method = 'BH',
            sort.by       = "M",
            resort.by     = "B",
            p.value       = 0.01,
            lfc           = 1,
            number        = Inf
        )
    tt$gene_name <- genes$SYMBOL[match(row.names(tt), genes$ENSEMBL)]
    write.csv(
        tt, 
        file = paste("./toptable_", coef, "_crh.csv", sep =""),
        row.names = TRUE
    )
    if (length(tt) > 0)    {
            plot_expression_tt(v, tt, crh.metadata, "group", "./", coef, 15)
    }

}
# ------------------- expression plot for any gene ----------------------

gene_list <- row.names(out)

# gene_list <- c("ENSMUSG00000036856", 
#                 "ENSMUSG00000032492", 
#                 "ENSMUSG00000044533", 
#                 "ENSMUSG00000038894", 
#                 "ENSMUSG00000020427", 
#                 "ENSMUSG00000028452",
#                 "ENSMUSG00000031906")
# gene_list <- rownames(v$E) %in% gene_list

pdf("./heatmap_custom_genelist_toptable_all_crh.pdf", height =  200, width = 16, onefile = TRUE)
    par(mar=c(5,6,4,2)+0.1); 
    x <- plot_expression(v, gene_list, crh.metadata, "group")
    print(x)
dev.off()

# de_gene_list_2 <- row.names(results[which(results[,1]!=0 | 
#                           results[,2]!=0 | 
#                           results[,3]!=0 | 
#                           results[,4]!=0 | 
#                           results[,5]!=0 | 
#                           results[,6]!=0 | 
#                           results[,7]!=0 | 
#                           results[,8]!=0 |
#                           results[,9]!=0 | 
#                           results[,10]!=0 | 
#                           results[,11]!=0 | 
#                           results[,12]!=0 | 
#                           results[,13]!=0 | 
#                           results[,14]!=0 | 
#                           results[,15]!=0 |
#                           results[,16]!=0 | 
#                           results[,17]!=0 | 
#                           results[,18]!=0 | 
#                           results[,19]!=0 | 
#                           results[,20]!=0
#                           ), ])

# pdf("./de_genelist_gene_level_heatmaps_all_p05.pdf", height =  70, width = 20, onefile = TRUE)
#     par(mar=c(5,6,4,2)+0.1); 
#     x <- plot_expression(v, de_gene_list, "expgroup")
#     print(x)
# dev.off()

# pdf("./de_genelist_gene_level_heatmaps_all_p1_dupcor.pdf", height =  50, width = 20, onefile = TRUE)
#     par(mar=c(5,6,4,2)+0.1); 
#     x <- plot_expression(v, de_gene_list_2, "expgroup")
#     print(x)
# dev.off()

# egr1 <- c("ENSMUSG00000038418")
# pdf("./egr1_gene_level_heatmaps_p05.pdf", height = 5, width = 20, onefile = TRUE)
#     par(mar=c(5,6,4,2)+0.1); 
#     x <- plot_expression(v, egr1, "expgroup")
#     print(x)
# dev.off()

# pappa2 <- c("ENSMUSG00000073530")
# pdf("./pappa_gene_level_heatmaps_p05.pdf", height = 5, width = 20, onefile = TRUE)
#     par(mar=c(5,6,4,2)+0.1); 
#     x <- plot_expression(v, pappa2, "expgroup")
#     print(x)
# dev.off()

# ------------------ Gene set enrichment (Not complete) -----------------------------

# MSigDB files
MM.c2 <- readRDS("../../MSigDB_v7.1/Mm.c2.all.v7.1.entrez.rds")
MM.c3 <- readRDS("../../MSigDB_v7.1/Mm.c3.all.v7.1.entrez.rds")
MM.c4 <- readRDS("../../MSigDB_v7.1/Mm.c4.all.v7.1.entrez.rds")
MM.c5 <- readRDS("../../MSigDB_v7.1/Mm.c5.all.v7.1.entrez.rds")
MM.c6 <- readRDS("../../MSigDB_v7.1/Mm.c6.all.v7.1.entrez.rds")
MM.c7 <- readRDS("../../MSigDB_v7.1/Mm.c7.all.v7.1.entrez.rds")
MM.h <- readRDS("../../MSigDB_v7.1/Mm.h.all.v7.1.entrez.rds")

# ------------------- Convert Ensemble to Entrez for DGE objects ----------

# retrieve the unique gene names in dge object
ensemble_to_entrez <- function(dge){
    symbolz           <- unique(dge$genes$genes)
    # use the org.Hs.eg.db to create a data frame of gene names and therir entrez IDs, add "NA" if not found.
    df_entrez         <- select(org.Mm.eg.db, symbolz, "ENTREZID","ENSEMBL")
    #Assigning the respective entrez ID to the dgeList object through matching.
    dge$genes$ENTREZ  <- df_entrez$ENTREZID[match(dge$genes$genes, df_entrez$ENSEMBL)]

    dge               <- dge[!is.na(dge$genes$ENTREZ),,keep.lib.sizes = FALSE]

    return(dge)
}

crh.dge <- ensemble_to_entrez(crh.dge)

# --------------------- create mroast objects ----------------------------

mroast.MM <- function(dge, efit, mm, design, mm_id, cont.matrix) {

    idx         <- limma::ids2indices(mm, unique(dge$genes$ENTREZ))
    # clean the lists to exclude "NA" values and gene sets which have less than 10 genes.
    idxclean    <- lapply(idx, function(x) x[!is.na(x)])
    idxlengths  <- sapply(idxclean, length)
    idxsub      <- idxclean[idxlengths > 10]

    for (contrast in colnames(cont.matrix)) {
        mroast_results <- mroast(
            dge,
            index    = idxsub,
            design   = design,
            contrast = cont.matrix[, contrast],
            nrot     = 99999
        )
        write.csv(mroast_results %>% tibble::rownames_to_column("Gene"), file = paste0("./GeneSetTesting/", contrast, "_", mm_id ,".csv"))
    }
}   

# Hallmark geneset
mroast.MM(crh.dge, efit, MM.h, crh.design, "hallmark", cont.matrix)
# Gene ontology 
mroast.MM(crh.dge, efit, MM.c5, crh.design, "GO", cont.matrix)
# KeGG
mroast.MM(crh.dge, efit, MM.c2, crh.design, "KEGG", cont.matrix)


pdf("../CRH/GeneSetTesting/UnxVsSham.hallmark.pdf", width = 10)
    hallmark_UvS <- read.csv("../CRH/GeneSetTesting/UnxvsSham_hallmark.csv", header = TRUE)
    p           <- plotTopGSEA (gseaout = hallmark_UvS, title = "GSEA Hallmark Unx Vs Sham", labCex = 2, textGap = 2)
    print(p)
dev.off()

pdf("../CRH/GeneSetTesting/UnxVsSham.GO.pdf", width = 100, heoght = 100)
    GO_UvS <- read.csv("../CRH/GeneSetTesting/UnxvsSham_GO.csv", header = TRUE)
    p      <- plotTopGSEA (gseaout = GO_UvS, title = "GSEA GO Unx Vs Sham", labCex = 2, textGap = 2)
    print(p)
dev.off()

# -------------------- Looking for Pappa2 regulators ---------------

# cd /fs02/ag36/Shani/Trev-Seq/
# less PAPPA2_regulators.R | tr ' ' '\n'| sort | uniq  > ../PAPPA2_regulators_edited.tab

pappa_regs      <- read.table("/fs02/ag36/Shani/Trev-Seq/PAPPA2_regulators_edited.tab")
pappa_reg_id_df <- subset(genes, genes$SYMBOL %in% tools::toTitleCase(tolower(pappa_regs$V1)))

pappa2          <- c("ENSMUSG00000073530")
egr1            <- c("ENSMUSG00000038418")

pdf("./heatmap_papparegs_genelist_toptable_all_crh.pdf", height =  200, width = 16, onefile = TRUE)
    par(mar=c(5,6,4,2)+0.1); 
    x <- plot_expression(v, c(pappa2, egr1, pappa_reg_id_df$ENSEMBL), crh.metadata, "group")
    print(x)
dev.off()

######################################################################
####################### helper functions #############################

# ----------------- Helper 01: expression plots for tt ----------------
plot_expression_tt <- function(v, tt, metadata, code, file.loc, coef, page_height = NULL) {
    avgexp                 <- subset(v$E, rownames(v$E) %in% rownames(tt))
    m.avgexp               <- melt(t(avgexp))
    colnames(m.avgexp)     <- c("ID","Gene","Expression")
    m.avgexp$timepoint     <- as.factor(metadata$timepoint[match(m.avgexp$ID,metadata$sampleName)])
    m.avgexp$treatment     <- as.factor(metadata$treatment[match(m.avgexp$ID,metadata$sampleName)])
    m.avgexp$label         <- as.factor(metadata$label[match(m.avgexp$ID,metadata$sampleName)])
    m.avgexp$group         <- as.factor(metadata$group[match(m.avgexp$ID,metadata$sampleName)])
    
    m.avgexp$common_name   <- tt$gene_name[match(m.avgexp$Gene, row.names(tt))]

    splitmeanE_gene       <- function(df, code) {
            s = split( df, df[, c(code, "Gene")]) #splits m.avgexp according to a given and gene
            sapply( s, function(x) mean(x$Expression) ) #calculates mean expression levels
        }

        u                     <- as.data.frame(splitmeanE_gene(m.avgexp, code))
        colnames(u)           <- c("Mean")
        m.avgexp$code.new     <- paste(m.avgexp$group, m.avgexp$Gene, sep=".")
        m.avgexp              <- cbind(u, 
                                        m.avgexp[match(row.names(u), m.avgexp$code.new),])
      	m.avgexp$Mean_transformed <- ifelse(m.avgexp$Mean<0.00,NA,m.avgexp$Mean)
		m.avgexp               <- plyr::ddply(m.avgexp, .(Gene), transform, rescl= scales::rescale(Mean_transformed))
        
        positive <- colorRampPalette(c("steelblue", "khaki1"))(50)                      
		negative <- colorRampPalette(c("white"))(1)
          g <- ggplot(
              m.avgexp,
              aes(
                  x    = group,
                  y    = factor(paste(m.avgexp$Gene, " (", m.avgexp$common_name,")" ,sep=""))
              )) +
              geom_tile(
                aes(
                  fill = rescl),
                colour = "black"
                ) +
              geom_text(
                aes(
                  label = round(Mean, 2)),
                  size  = 5
                )+
              scale_fill_gradientn(
				colours=c(positive), 
				na.value = "white",
				limits = c(0, 1), 
                guide = "none"
                ) + 
               theme_minimal() +
               theme(
                  strip.text          = element_text(size=24, face = "bold"), #base_size didin't work
                  plot.title          = element_text(size = 28, vjust=1, hjust=1, face = "bold"),
                  panel.grid.major    = element_line(size=0.1),
                  panel.grid.minor    = element_line(size=0.5),
                  panel.spacing       = unit(1.5, 'lines'),
                  axis.title.y        = element_text(angle=90, vjust=1.25, size=1, face="bold"),
                  axis.title.x        = element_text(vjust=1, size=1, face="bold"),
                  axis.text.y         = element_text(size=14, face="italic"),
                  axis.text.x         = element_text(angle = 90, hjust = 1, size=24)
              )+
              scale_x_discrete(name="",expand = c(0, 0)) +
              scale_y_discrete(name="",expand = c(0, 0)) +  
              facet_grid(~timepoint, scales='free_x', space="free_x") + 
              ggtitle(coef)
		
        if (is.null(page_height)) {page_height <- 10}
        pdf(paste(file.loc, coef, "_gene_level_heatmap_crh.pdf", sep =""), height =  page_height, width = 20)
            par(mar=c(5,6,4,2)+0.1); 
            print(g)
        dev.off()
}

# ----------------- Helper 02: expression plots for custom set of genes ----------------
plot_expression <- function(v, gene_subset, metadata, code) {
    avgexp                 <- subset(v$E, rownames(v$E) %in% gene_subset)
    m.avgexp               <- melt(t(avgexp))
    colnames(m.avgexp)     <- c("ID","Gene","Expression")
    m.avgexp$timepoint     <- as.factor(metadata$timepoint[match(m.avgexp$ID,metadata$sampleName)])
    m.avgexp$treatment     <- as.factor(metadata$treatment[match(m.avgexp$ID,metadata$sampleName)])
    m.avgexp$label         <- as.factor(metadata$label[match(m.avgexp$ID,metadata$sampleName)])
    m.avgexp$group         <- as.factor(metadata$group[match(m.avgexp$ID,metadata$sampleName)])

    m.avgexp$common_name   <- genes$SYMBOL[match(m.avgexp$Gene, genes$ENSEMBL)]

    splitmeanE_gene       <- function(df, code) {
            s = split( df, df[, c(code, "Gene")]) #splits m.avgexp according to a given and gene
            sapply( s, function(x) mean(x$Expression) ) #calculates mean expression levels
        }

        u                     <- as.data.frame(splitmeanE_gene(m.avgexp, code))
        colnames(u)           <- c("Mean")
        m.avgexp$code.new     <- paste(m.avgexp$group, m.avgexp$Gene, sep=".")
        m.avgexp              <- cbind(u, 
                                        m.avgexp[match(row.names(u), m.avgexp$code.new),])
      	m.avgexp$Mean_transformed <- ifelse(m.avgexp$Mean<0.00,NA,m.avgexp$Mean)
		m.avgexp               <- plyr::ddply(m.avgexp, .(Gene), transform, rescl= scales::rescale(Mean_transformed))
        
        positive <- colorRampPalette(c("steelblue", "khaki1"))(50)                      
		negative <- colorRampPalette(c("white"))(1)
          g <- ggplot(
              m.avgexp,
              aes(
                  x    = group,
                  y    = factor(paste(m.avgexp$Gene, " (", m.avgexp$common_name,")" ,sep=""))
              )) +
              geom_tile(
                aes(
                  fill = rescl),
                colour = "black"
                ) +
              geom_text(
                aes(
                  label = round(Mean, 2)),
                  size  = 4
                )+
              scale_fill_gradientn(
				colours=c(positive), 
				na.value = "white",
				limits = c(0, 1), 
                guide = "none"
                ) + 
               theme_minimal() +
               theme(
                  strip.text          = element_text(size=24, face = "bold"), #base_size didin't work
                  plot.title          = element_text(size = 28, vjust=1, hjust=1, face = "bold"),
                  panel.grid.major    = element_line(size=0.1),
                  panel.grid.minor    = element_line(size=0.5),
                  panel.spacing       = unit(1.5, 'lines'),
                  axis.title.y        = element_text(angle=90, vjust=1.25, size=1, face="bold"),
                  axis.title.x        = element_text(vjust=1, size=1, face="bold"),
                  axis.text.y         = element_text(size=14, face="italic"),
                  axis.text.x         = element_text(angle = 90, hjust = 1, size=24)
              )+
              scale_x_discrete(name="",expand = c(0, 0)) +
              scale_y_discrete(name="",expand = c(0, 0)) +  
              facet_grid(~timepoint, scales='free_x', space="free_x")
    g
}

# ------------------ INTERACTIVE VOLCANO PLOT --------------------

volcano_plot_interactive <- function(efit, dge, main, html) {
  cols <- if (length(unique(decideTests(efit))) > 1) {
    c("orange", "black", "orange")
  } else {
    c("#00bfff", "#858585", "#00bfff")
  }

  colourCount  <- length(unique(dge$sample$group))
  getPalette   <- colorRampPalette(brewer.pal(6, "Dark2"))
  col_groups   <- getPalette(colourCount)

  glXYPlot(
    efit$coefficient,
    -log10(efit$p.value),
    main        = main,
    xlab        = "log2-fold-change",
    ylab        = "-log10(p-value)",
    counts      = cpm(dge$counts),
    side.main   = "genes",
    anno        = dge$genes,
    group       = as.factor(as.character(dge$sample$group)),
    status      = decideTests(efit),
    cols        = cols,
    sample.cols = col_groups[as.factor(as.character(dge$sample$group))],
    html        = html,
    folder      = "./glimma_plots",
    launch      = FALSE
  )
}

# ------------------INTEREACTIVE MD PLOT -----------------

md_plot_interactive <- function(efit, dge, coef.id, html) {
#   cols <- if (length(unique(decideTests(efit))) > 1) {
#     c("orange", "black", "orange")
#   } else {
#     c("#00bfff", "#858585", "#00bfff")
#   }
    colourCount  <- length(unique(dge$sample$group))
    getPalette   <- colorRampPalette(brewer.pal(6, "Dark2"))
    col_groups   <- getPalette(colourCount)
  glMDPlot(efit, 
           coef        = coef.id, 
           status      = decideTests(efit),
           main        = colnames(efit)[coef.id],
           counts      = cpm(dge$counts),
           anno        = dge$genes,
           side.main   = "ENSEMBL",
           group       = dge$sample$group,
           sample.cols = col_groups[as.factor(as.character(dge$sample$group))],
           html        = html,
           folder      = "./glimma_plots",
           launch      = FALSE
           )
}

# --------------- Plotting GSEA results -----------------

plotTopGSEA <- function(gseaout, topN = NULL, title,
                         gapProp = 0.5, textGap = 5, textIndent = 1, labCex = 3){
  require(ggplot2)

  gseaout$PropDown <- round(gseaout$PropDown, 3)
  gseaout$PropUp   <- round(gseaout$PropUp, 3)
  if (is.null(topN)) {
    topN <- NROW(gseaout)
  }
  
  GOtopT <- gseaout[1:topN, ]
  offset <- ifelse((sum(log(GOtopT$FDR)==Inf)+sum(log(GOtopT$FDR)==-Inf))>0, 10^(-10), 100)
  labels <- paste0(GOtopT$Gene)[topN:1]
  d <- data.frame(x = c(topN:1, topN:1),
                  Direction = c(rep("down",topN), rep("up",topN)),
                  logP = c(ifelse(GOtopT$Direction=="Down", log(GOtopT$PValue)*offset, 0), ifelse(GOtopT$Direction=="Up", -log(GOtopT$PValue)*offset, 0)), #c(log(GOtopT$PropDown + offset), -log(GOtopT$PropUp + offset)), 
                  n = c(GOtopT$NGenes, GOtopT$NGenes),
                  nDirection = c(GOtopT$PropDown, GOtopT$PropUp))
  #
  axisGap   <- (max(d$logP) - min(d$logP))*gapProp
  VlineGap  <- axisGap*0.25
  topLabGap <- axisGap*0.5
  Ymin <- min(c(-topLabGap, min(d$logP), log(0.05 + offset))) - axisGap
  Ymax <- max(c(topLabGap, max(d$logP), -log(0.05 + offset))) + textGap  
  #
  p <- ggplot(d, aes(x=x, y=logP, fill=Direction)) +
    geom_bar(data = subset(d, Direction == "down"),
             aes(y = logP), stat = "identity", position = "dodge", color = "lightblue", fill="lightblue") +
    geom_bar(data = subset(d, Direction == "up"),
             aes(y = logP), stat = "identity", position = "dodge", color = "lightpink", fill="lightpink") +
    geom_hline(yintercept = 0,colour = "grey90")+
    coord_flip() +
    theme_classic() +
    scale_x_continuous(breaks=1:topN, labels=labels)+
    ylim(c(Ymin, Ymax)) +
    labs(title= title) +
    theme(axis.title.y=element_blank(),
          axis.title.x = element_blank(),
          plot.title = element_text(hjust = 0.5),
          legend.position = "none") +
    annotate("text", x = topN+1, y = Ymin+VlineGap-textIndent,
             label="n", color = "black", hjust = "right", size=labCex) +
    annotate("text", x = topN+1, y = Ymin+(2*VlineGap)-textIndent,
             label="Down", color = "lightblue", hjust = "right", size=labCex) +
    annotate("text", x = topN+1, y = Ymin+(3*VlineGap)-textIndent,
             label="Up", color = "lightpink", hjust = "right", size=labCex) +
    annotate("text", x = topN+1, y = -topLabGap,
             label="log(P.Down)", color = "lightblue", size=labCex) +
    annotate("text", x = topN+1, y = topLabGap,
             label="-log(P.Up)", color = "lightpink", size=labCex) +
    geom_text(data = subset(d, Direction == "down"),
              aes(y = Ymin+VlineGap-textIndent, label=n),
              stat = "identity", color = "black", hjust = "right", size=labCex) +
    geom_text(data = subset(d, Direction == "down"),
              aes(y = Ymin+(2*VlineGap)-textIndent, label=nDirection),
              stat = "identity", color = "lightblue", hjust = "right", size=labCex) +
    geom_text(data = subset(d, Direction == "up"),
              aes(y = Ymin+(3*VlineGap)-textIndent, label=nDirection),
              stat = "identity", color = "lightpink", hjust = "right", size=labCex)  +
    geom_hline(yintercept = Ymin+VlineGap, colour = "gray") +
    geom_hline(yintercept = Ymin+(2*VlineGap), colour = "gray") +
    geom_hline(yintercept = Ymin+(3*VlineGap), colour = "gray") +
    geom_vline(xintercept = topN+0.5, colour = "gray") +
    #FDR 0.1 lines:
    geom_hline(yintercept = log(0.1)*offset, colour = "black", linetype = "dashed") +
    geom_hline(yintercept = -log(0.1)*offset, colour = "black", linetype = "dashed")
  
    return(p)    
}