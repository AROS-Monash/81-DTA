library(msigdbr)
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
library(tidyverse)
library(magrittr)
library(viridis)


hm <- msigdbr("Mus musculus", category = "H") %>% 
                dplyr::distinct(gs_name, entrez_gene, .keep_all = TRUE) %>% 
                as.data.frame()


hm <- msigdbr("Mus musculus", category = "H")  %>% 
  #left_join(entrez_gene) %>%
  dplyr::filter(!is.na(ensembl_gene)) %>%
  distinct(gs_name, ensembl_gene, .keep_all = TRUE)
hmByGene <- hm %>%
  split(f = .$ensembl_gene) %>%
  lapply(extract2, "gs_name")
hmByID <- hm %>%
  split(f = .$gs_name) %>%
  lapply(extract2, "ensembl_gene")

#there is VAX and IMMUNESIGDB, we picked IMMUNESIGDB only
imm <- msigdbr("Mus musculus", category = "C7", subcategory = "IMMUNESIGDB")  %>% 
  #left_join(entrezGenes) %>%
  dplyr::filter(!is.na(ensembl_gene)) %>%
  distinct(gs_name, ensembl_gene, .keep_all = TRUE)
immByGene <- imm %>%
  split(f = .$ensembl_gene) %>%
  lapply(extract2, "gs_name")
immByID <- imm %>%
  split(f = .$gs_name) %>%
  lapply(extract2, "ensembl_gene")

kg <- msigdbr("Mus musculus", category = "C2", subcategory = "CP:KEGG")  %>% 
  #left_join(entrezGenes) %>%
  dplyr::filter(!is.na(ensembl_gene)) %>%
  distinct(gs_name, ensembl_gene, .keep_all = TRUE)
kgByGene <- kg %>%
  split(f = .$ensembl_gene) %>%
  lapply(extract2, "gs_name")
kgByID <- kg %>%
  split(f = .$gs_name) %>%
  lapply(extract2, "ensembl_gene")

TF <- msigdbr("Mus musculus", category = "C3")  %>% 
  dplyr::filter(gs_subcat==c("TFT:GTRD","TFT:TFT_Legacy")) %>%
  #left_join(entrezGenes) %>%
  dplyr::filter(!is.na(ensembl_gene)) %>%
  distinct(gs_name, ensembl_gene, .keep_all = TRUE)
TFByGene <- TF %>%
  split(f = .$ensembl_gene) %>%
  lapply(extract2, "gs_name")
TFByID <- TF %>%
  split(f = .$gs_name) %>%
  lapply(extract2, "ensembl_gene")

miR <- msigdbr("Mus musculus", category = "C3")  %>% 
  dplyr::filter(gs_subcat==c("MIR:MIR_Legacy","MIR:MIRDB")) %>%
  #left_join(entrezGenes) %>%
  dplyr::filter(!is.na(ensembl_gene)) %>%
  distinct(gs_name, ensembl_gene, .keep_all = TRUE)
miRByGene <- miR%>%
  split(f = .$ensembl_gene) %>%
  lapply(extract2, "gs_name")
miRByID <- miR %>%
  split(f = .$gs_name) %>%
  lapply(extract2, "ensembl_gene")

######## Adelaide E15.5 ################

setwd(/home/samarasi/ag36_scratch/Shani/Trev-Seq/Adelaide_approach/E15)

samples <- readRDS(here("data", "samples.rds"))
dgeList <- readRDS(here("data", "dgeList.rds"))
topTables <- readRDS(here("data", "topTables.rds"))
v <- readRDS(here("data", "v.rds"))

design <- v$design

contr.matrix <- makeContrasts(
   ExpLvsExpR = Exp.L-Exp.R, 
   ExpLvsCtrl= Exp.L-Ctrl, 
   ExpRvsCtrl = Exp.R-Ctrl,
   ExpvsCtrl = (Exp.L+Exp.R)-Ctrl, 
   levels = colnames(design))
contr.matrix

entrezGenes <- dgeList$genes %>%
  dplyr::filter(!is.na(entrezid)) %>%
  tidyr::unnest(entrezid) %>%
  dplyr::rename(entrez_gene = entrezid) %>%
  mutate(entrez_gene=as.integer(as.character(entrez_gene)))

deTable <- names(topTables) %>% 
  sapply(function(x){ topTables[[x]] %>%
      dplyr::filter() %>%
  mutate(entrezid = dgeList$genes$entrezid[gene_id])}, simplify = FALSE)

hm <- msigdbr("Mus musculus", category = "H")  %>% 
  left_join(entrezGenes) %>%
  dplyr::filter(!is.na(gene_id)) %>%
  distinct(gs_name, gene_id, .keep_all = TRUE)
hmByGene <- hm %>%
  split(f = .$gene_id) %>%
  lapply(extract2, "gs_name")
hmByID <- hm %>%
  split(f = .$gs_name) %>%
  lapply(extract2, "gene_id")

#there is VAX and IMMUNESIGDB, we picked IMMUNESIGDB only
imm <- msigdbr("Mus musculus", category = "C7", subcategory = "IMMUNESIGDB")  %>% 
  left_join(entrezGenes) %>%
  dplyr::filter(!is.na(gene_id)) %>%
  distinct(gs_name, gene_id, .keep_all = TRUE)
immByGene <- imm %>%
  split(f = .$gene_id) %>%
  lapply(extract2, "gs_name")
immByID <- imm %>%
  split(f = .$gs_name) %>%
  lapply(extract2, "gene_id")

kg <- msigdbr("Mus musculus", category = "C2", subcategory = "CP:KEGG")  %>% 
  left_join(entrezGenes) %>%
  dplyr::filter(!is.na(gene_id)) %>%
  distinct(gs_name, gene_id, .keep_all = TRUE)
kgByGene <- kg %>%
  split(f = .$gene_id) %>%
  lapply(extract2, "gs_name")
kgByID <- kg %>%
  split(f = .$gs_name) %>%
  lapply(extract2, "gene_id")

TF <- msigdbr("Mus musculus", category = "C3")  %>% 
  dplyr::filter(gs_subcat==c("TFT:GTRD","TFT:TFT_Legacy")) %>%
  left_join(entrezGenes) %>%
  dplyr::filter(!is.na(gene_id)) %>%
  distinct(gs_name, gene_id, .keep_all = TRUE)
TFByGene <- TF %>%
  split(f = .$gene_id) %>%
  lapply(extract2, "gs_name")
TFByID <- TF %>%
  split(f = .$gs_name) %>%
  lapply(extract2, "gene_id")

miR <- msigdbr("Mus musculus", category = "C3")  %>% 
  dplyr::filter(gs_subcat==c("MIR:MIR_Legacy","MIR:MIRDB")) %>%
  left_join(entrezGenes) %>%
  dplyr::filter(!is.na(gene_id)) %>%
  distinct(gs_name, gene_id, .keep_all = TRUE)
miRByGene <- miR%>%
  split(f = .$gene_id) %>%
  lapply(extract2, "gs_name")
miRByID <- miR %>%
  split(f = .$gs_name) %>%
  lapply(extract2, "gene_id")

gsSizes <- names(deTable) %>% sapply(function(x){bind_rows(hm, kg, TF, miR, imm) %>%
  dplyr::select(gs_name, gene_symbol, gene_id, gs_id, gs_pmid, gs_geoid) %>% 
  chop(c(gene_symbol, gene_id)) %>%
  mutate(
    gs_size = vapply(gene_symbol, length, integer(1)),
    de_id = lapply(
      X = gene_id, 
      FUN = intersect, 
      y = dplyr::filter(deTable[[x]], DE)$gene_id
      ),
    de_size = vapply(de_id, length, integer(1))
  )}, simplify=FALSE)


hmFry <- names(topTables) %>% 
  sapply(function(x){ v %>%
  fry(index = hmByID,
    design = design,
    contrast = contr.matrix[,x],
    sort = "directional") %>%
  rownames_to_column("gs_name") %>%
  as_tibble()%>%
  left_join(gsSizes[[x]]) %>%
  mutate(gs_name = str_remove(gs_name, "HALLMARK_"),
         gs_name = gsub("_", " ",gs_name)) }, simplify = FALSE)

saveRDS(hmFry, file="./data/hmFry.rds")

#this step is select top 10 Hallmark from each group
# PLEASE CHANGE FDR cutoff HERE!!!!!
hmgroup<-names(hmFry)%>% 
  sapply(function(x){hmFry[[x]] %>% dplyr::slice(1:30) %>% .$gs_name}, simplify = FALSE)
# combine all and create a unique group of hm that are sig
hmgroupTOP<-c(hmgroup[[1]],hmgroup[[2]],hmgroup[[3]],hmgroup[[4]]) %>% unique()

hm_top30 <-names(hmFry) %>% 
  sapply(function(x){hmFry[[x]] %>% 
  dplyr::filter(gs_name %in% hmgroupTOP) %>%
  mutate(gene_ratio=de_size/NGenes) %>%
  dplyr::select("gs_name","NGenes","Direction","PValue", "FDR","de_size","gene_ratio") %>%
  mutate(Group=x)}, simplify = FALSE)

# Arrange the hallmark genesets based on level 
hmgroupTOPlvl <- hm_top30$ExpLvsExpR$gs_name  %>% factor(levels=.)

hmALL<-rbind(hm_top30$ExpLvsExpR, hm_top30$ExpLvsCtrl,hm_top30$ExpRvsCtrl) %>%
  mutate(gs_name=factor(gs_name,levels=levels(fct_rev(hmgroupTOPlvl)))) %>%
  dplyr::select(Hallmark="gs_name", everything()) %>%
  mutate(Group=factor(Group, 
                       levels=c("ExpLvsExpR", "ExpLvsCtrl","ExpRvsCtrl"))) 

saveRDS(hmALL,file="./data/hmALL.rds")

# hmplotting
hmALL_plot <- hmALL %>% ggplot(aes(x= Group, y=Hallmark,  color=gene_ratio, shape=Direction)) + 
  geom_point(aes(size=-log10(PValue))) +  
  scale_color_viridis(option="B", direction=-1, begin=0.4, end= 0.8)+
  scale_size(range=c(1,7)) + #adjust the size of the circle
  theme_classic() +
  theme(axis.text.x=element_text(angle = 45,vjust = 1, hjust=1),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        panel.grid.major = element_line(colour="gray90"))+
    labs(size= "-log10(PValue)", colour="Gene ratio")


pdf("./hmAll_plot.pdf", height = 10)
hmALL_plot
dev.off()








##### E15.5 ###########

setwd("/fs03/ag36/Shani/Trev-Seq/featurecounts_edgeR_Mar_2023/E15")

# save the v.15, dge.15 and efit.15
dge.15    <- readRDS("./dge15.final.rds")
v.15      <- readRDS("./voom15.final.rds")
efit.15   <- readRDS("./efit15.final.rds")
design.15 <- readRDS("./design15.final.rds")

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

gsSizes <- colnames(efit.15) %>% sapply(function(x){bind_rows(hm, kg, TF, miR, imm) %>%
  dplyr::select(gs_name, ensembl_gene, gs_id, gs_pmid, gs_geoid, gene_symbol) %>% 
  chop(c(ensembl_gene, gene_symbol)) %>%
  mutate(
    gs_size = vapply(ensembl_gene, length, integer(1)),
    de_id = lapply(
      X = ensembl_gene, 
      FUN = intersect, 
      y = row.names(efit.15$genes)
      ),
    de_size = vapply(de_id, length, integer(1))
  )}, simplify=FALSE)

hmFry.15 <- colnames(efit.15) %>% 
  sapply(function(x){ v.15 %>%
  fry(index = hmByID,
    design = design.15,
    contrast = cont.matrix.15[,x],
    sort = "directional") %>%
  rownames_to_column("gs_name") %>%
  as_tibble()%>%
  left_join(gsSizes[[x]]) %>%
  mutate(gs_name = str_remove(gs_name, "HALLMARK_"),
         gs_name = gsub("_", " ",gs_name)) }, simplify = FALSE)

saveRDS(hmFry.15, file="./GeneSetTesting/hmFry.15.rds")

#this step is select top 10 Hallmark from each group
# PLEASE CHANGE FDR cutoff HERE!!!!!
hmgroup<-names(hmFry.15)%>% 
  sapply(function(x){hmFry.15[[x]] %>% filter(P.Value < 0.1) %>% dplyr::slice(1:30) %>% .$gs_name}, simplify = FALSE)
# combine all and create a unique group of hm that are sig
hmgroupTOP<-c(hmgroup[[1]],hmgroup[[2]],hmgroup[[3]],hmgroup[[4]], hmgroup[[5]], hmgroup[[6]]) %>% unique() 

hm_top30 <-names(hmFry.15) %>% 
  sapply(function(x){hmFry.15[[x]] %>% 
  dplyr::filter(gs_name %in% hmgroupTOP) %>%
  mutate(gene_ratio=de_size/NGenes) %>%
  dplyr::select("gs_name","NGenes","Direction","PValue", "FDR","de_size","gene_ratio") %>%
  mutate(Group=x)}, simplify = FALSE)

# Arrange the hallmark genesets based on level 
hmgroupTOPlvl <- hm_top30$E_L_15vsRest$gs_name  #%>% factor(levels=.)

hmALL<-rbind(hm_top30$E_L_15vsRest, hm_top30$E_R_15vsRest,hm_top30$E_L_15vsE_R_15) %>%
  mutate(gs_name=factor(gs_name,levels=levels(fct_rev(hmgroupTOPlvl)))) %>%
  dplyr::select(Hallmark="gs_name", everything()) %>%
  mutate(Group=factor(Group, 
                       levels=c("E_L_15vsRest", "E_R_15vsRest","E_L_15vsE_R_15"))) 

saveRDS(hmALL,file="./GeneSetTesting/hmALL_15.rds")

# hmplotting
hmALL_plot <- hmALL %>% ggplot(aes(x= Group, y=Hallmark,  color=gene_ratio, shape=Direction)) + 
  geom_point(aes(size=-log10(PValue))) +  
  scale_color_viridis(option="B", direction=-1, begin=0.4, end= 0.8)+
  scale_size(range=c(1,7)) + #adjust the size of the circle
  theme_classic() +
  theme(axis.text.x=element_text(angle = 45,vjust = 1, hjust=1),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        panel.grid.major = element_line(colour="gray90"))+
    labs(size= "-log10(PValue)", colour="Gene ratio")

pdf("msigdbr_hallmark_15.pdf")
    hmALL_plot
dev.off()
########## extra

# for (coef in colnames(efit.15)) {
#     tt <- topTable(
#         efit.15,
#         coef = coef,
#         adjust.method = 'BH',
#         sort.by = "M",
#         resort.by = "B",
#         p.value = 0.1,
#         lfc = log2(1),
#         number = Inf
#     )
#     if (nrow(tt) > 0)    {
#         enricher(gene = row.names(tt), TERM2GENE = hm)
#     }
# }

# # E15.5

# deTable <- names(tt) %>% 
#   sapply(function(x){ tt[[x]] %>%
#       dplyr::filter(), , simplify = FALSE)

gsSizes <- colnames(efit.15) %>% sapply(function(x){bind_rows(hm, kg, TF, miR, imm) %>%
  dplyr::select(gs_name, ensembl_gene, gs_id, gs_pmid, gs_geoid, gene_symbol) %>% 
  chop(c(ensembl_gene, gene_symbol)) %>%
  mutate(
    gs_size = vapply(ensembl_gene, length, integer(1)),
    de_id = lapply(
      X = ensembl_gene, 
      FUN = intersect, 
      y = row.names(efit.15$genes)
      ),
    de_size = vapply(de_id, length, integer(1))
  )}, simplify=FALSE)