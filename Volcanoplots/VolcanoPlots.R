#
#Script for generating VolcanoPlots after DE and highlighting common genes among mice
#
#
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggrepel)
#
#
#
#
#
#Create the pdf
pdf("VolcanoPlots.pdf")
#
#We first load the .csv files for all three acute datasets (HIV+ high vs. HIV-)
MAST_297 <- read.table('HIVhighpos-vs.-HIVneg-297_count-MAST-markers.csv', header = TRUE, sep = ',', stringsAsFactors = FALSE)
MAST_618 <- read.table('HIVhighpos-vs.-HIVneg-618_count-MAST-markers.csv', header = TRUE, sep = ',', stringsAsFactors = FALSE)
MAST_620 <- read.table('HIVhighpos-vs.-HIVneg-620_count-MAST-markers.csv', header = TRUE, sep = ',', stringsAsFactors = FALSE)

MAST_297<-rename(MAST_297, "Gene" = "X")
MAST_618<-rename(MAST_618, "Gene" = "X")
MAST_620<-rename(MAST_620, "Gene" = "X")


#We make a copy with only those that are fdr<0.05
MAST_297$threshold <- ifelse(MAST_297$fdr < 0.05, "FDR < 0.05", "Not Sig")
MAST_618$threshold <- ifelse(MAST_618$fdr < 0.05, "FDR < 0.05", "Not Sig")
MAST_620$threshold <- ifelse(MAST_620$fdr < 0.05, "FDR < 0.05", "Not Sig")

#
MAST_297_sig<-subset.data.frame(MAST_297, threshold == "FDR < 0.05")
MAST_618_sig<-subset.data.frame(MAST_618, threshold == "FDR < 0.05")
MAST_620_sig<-subset.data.frame(MAST_620, threshold == "FDR < 0.05")

#We make copies and only take the genes columns

MAST_297_G<-subset.data.frame(MAST_297_sig,select = "Gene")
MAST_618_G<-subset.data.frame(MAST_618_sig,select = "Gene")
MAST_620_G<-subset.data.frame(MAST_620_sig,select = "Gene")


# Join Acute Genes
joint.bcs <- inner_join(MAST_297_G, MAST_618_G)
common_genes <- inner_join(joint.bcs, MAST_620_G)

#Make it into a vector
common_genes <- common_genes$Gene


#Volcano plots using ggplot (HIV+ high vs. HIV-) -A1
MAST_297$Commonality <- ifelse(MAST_297$Gene %in% common_genes, "FDR < 0.05 in all three", "Not")

MAST_297$vplot <- case_when( MAST_297$Commonality == "FDR < 0.05 in all three" ~ "FDR < 0.05 in all three",
                             MAST_297$threshold == "FDR < 0.05" ~ "FDR < 0.05",.default = "Not Sig")


ggplot(MAST_297, aes(x = avg_log2FC, y = -log10(fdr))) +
  geom_point(aes(color = vplot)) +
  scale_color_manual(values = c("red", "blue","lightgray")) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom") +
  ggtitle("A1 15dpi HIV+ High vs HIV- Volcano Plot") +
  coord_cartesian(xlim = c(-4,4))+
  geom_text_repel(data = subset(MAST_297, vplot == "FDR < 0.05 in all three"),
    aes(label = Gene),size = 3,box.padding = unit(0.1, "lines"),point.padding = unit(0.1, "lines")) 


#Volcano plots using ggplot (HIV+ high vs. HIV-) -A2

MAST_618$Commonality <- ifelse(MAST_618$Gene %in% common_genes, "FDR < 0.05 in all three", "Not")


MAST_618$vplot <- case_when( MAST_618$Commonality == "FDR < 0.05 in all three" ~ "FDR < 0.05 in all three",
                             MAST_618$threshold == "FDR < 0.05" ~ "FDR < 0.05",.default = "Not Sig")



ggplot(MAST_618, aes(x = avg_log2FC, y = -log10(fdr))) +
  geom_point(aes(color = vplot)) +
  scale_color_manual(values = c("red", "blue","lightgray")) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom") +
  ggtitle("A2 15dpi HIV+ High vs HIV- Volcano Plot") +
  coord_cartesian(xlim = c(-5,5))+
  geom_text_repel(
    data = subset(MAST_618, vplot == "FDR < 0.05 in all three"),
    aes(label = Gene),
    size = 3,
    box.padding = unit(0.1, "lines"),
    point.padding = unit(0.1, "lines")
  ) 

#Volcano plots using ggplot (HIV+ high vs. HIV-) - A3

MAST_620$Commonality <- ifelse(MAST_620$Gene %in% common_genes, "FDR < 0.05 in all three", "Not")

MAST_620$vplot <- case_when( MAST_620$Commonality == "FDR < 0.05 in all three" ~ "FDR < 0.05 in all three",
                             MAST_620$threshold == "FDR < 0.05" ~ "FDR < 0.05",.default = "Not Sig")



ggplot(MAST_620, aes(x = avg_log2FC, y = -log10(fdr))) +
  geom_point(aes(color = vplot)) +
  scale_color_manual(values = c("red", "blue","lightgray")) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom") +
  ggtitle("A3 15dpi HIV+ High vs HIV- Volcano Plot") +
  coord_cartesian(xlim = c(-5,5))+
  geom_text_repel(
    data = subset(MAST_620, vplot == "FDR < 0.05 in all three"),
    aes(label = Gene),
    size = 3,
    box.padding = unit(0.1, "lines"),
    point.padding = unit(0.1, "lines")
  ) 

#####
#### Volcano plot for GFP vs. dsRed (T1/T2 , A2 ,A3)

MAST_368_GFP <- read.table('dsRed-vs.-GFP-368_count-MAST-markers.csv', header = TRUE, sep = ',', stringsAsFactors = FALSE)
MAST_618_GFP <- read.table('dsRed-vs.-GFP-618_count-MAST-markers.csv', header = TRUE, sep = ',', stringsAsFactors = FALSE)
MAST_620_GFP <- read.table('dsRed-vs.-GFP-620_count-MAST-markers.csv', header = TRUE, sep = ',', stringsAsFactors = FALSE)


MAST_368_GFP<-rename(MAST_368_GFP, "Gene" = "X")
MAST_618_GFP<-rename(MAST_618_GFP, "Gene" = "X")
MAST_620_GFP<-rename(MAST_620_GFP, "Gene" = "X")


#We make a copy with only those that are fdr<0.05
MAST_368_GFP$threshold <- ifelse(MAST_368_GFP$fdr < 0.05, "FDR < 0.05", "Not Sig")
MAST_618_GFP$threshold <- ifelse(MAST_618_GFP$fdr < 0.05, "FDR < 0.05", "Not Sig")
MAST_620_GFP$threshold <- ifelse(MAST_620_GFP$fdr < 0.05, "FDR < 0.05", "Not Sig")

#
MAST_368_GFP_sig<-subset.data.frame(MAST_368_GFP, threshold == "FDR < 0.05")
MAST_618_GFP_sig<-subset.data.frame(MAST_618_GFP, threshold == "FDR < 0.05")
MAST_620_GFP_sig<-subset.data.frame(MAST_620_GFP, threshold == "FDR < 0.05")

#We make copies and only take the genes columns

MAST_368_GFP_G<-subset.data.frame(MAST_368_GFP_sig,select = "Gene")
MAST_618_GFP_G<-subset.data.frame(MAST_618_GFP_sig,select = "Gene")
MAST_620_GFP_G<-subset.data.frame(MAST_620_GFP_sig,select = "Gene")



# Join Acute Genes
joint.bcs.GFP <- inner_join(MAST_368_GFP_G, MAST_618_GFP_G)
common_genes.GFP <- inner_join(joint.bcs.GFP, MAST_620_GFP_G)

#Make it into a vector
common_genes.GFP <- common_genes.GFP$Gene


#Volcano plot using ggplot (T1/T2)

MAST_368_GFP$Commonality<-ifelse(MAST_368_GFP$Gene %in% common_genes.GFP, "FDR < 0.05 in all three", "Not")

MAST_368_GFP$vplot <- case_when(MAST_368_GFP$Commonality == "FDR < 0.05 in all three" ~ "FDR < 0.05 in all three",
                                MAST_368_GFP$threshold == "FDR < 0.05" ~ "FDR < 0.05",.default = "Not Sig")


ggplot()+
  geom_point(data = MAST_368_GFP[MAST_368_GFP$vplot == "Not Sig",], aes(x = avg_log2FC, y = -log10(fdr),color=vplot))+
  geom_point(data = MAST_368_GFP[MAST_368_GFP$vplot == "FDR < 0.05",], aes(x = avg_log2FC, y = -log10(fdr),color=vplot))+
  geom_point(data = MAST_368_GFP[MAST_368_GFP$vplot == "FDR < 0.05 in all three",], aes(x = avg_log2FC, y = -log10(fdr),color=vplot))+
  scale_color_manual(values = c("red", "darkgreen", "lightgray")) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom") +
  ggtitle("T1/T2 10dpt/29dpt GFP vs dsRed Volcano Plot") +
  coord_cartesian(xlim = c(-9,9))+
  geom_text_repel(data = subset(MAST_368_GFP, vplot == "FDR < 0.05 in all three"),
                  aes(label = Gene,x = avg_log2FC, y = -log10(fdr)),
                  size = 3,max.overlaps = 20
  ) 



#Volcano plot using ggplot -A2

MAST_618_GFP$Commonality<-ifelse(MAST_618_GFP$Gene %in% common_genes.GFP, "FDR < 0.05 in all three", "Not")

MAST_618_GFP$vplot <- case_when(MAST_618_GFP$Commonality == "FDR < 0.05 in all three" ~ "FDR < 0.05 in all three",
                                MAST_618_GFP$threshold == "FDR < 0.05" ~ "FDR < 0.05",.default = "Not Sig")


ggplot()+
  geom_point(data = MAST_618_GFP[MAST_618_GFP$vplot == "Not Sig",], aes(x = avg_log2FC, y = -log10(fdr),color=vplot))+
  geom_point(data = MAST_618_GFP[MAST_618_GFP$vplot == "FDR < 0.05",], aes(x = avg_log2FC, y = -log10(fdr),color=vplot))+
  geom_point(data = MAST_618_GFP[MAST_618_GFP$vplot == "FDR < 0.05 in all three",], aes(x = avg_log2FC, y = -log10(fdr),color=vplot))+
  scale_color_manual(values = c("red", "darkgreen", "lightgray")) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom") +
  ggtitle("A2 15dpi GFP vs dsRed Volcano Plot") +
  geom_text_repel(data = subset(MAST_618_GFP, vplot == "FDR < 0.05 in all three"),
                  aes(label = Gene,x = avg_log2FC, y = -log10(fdr)),
                  size = 3,max.overlaps = 20
  ) 


#Volcano plot using ggplot - A3

MAST_620_GFP$Commonality<-ifelse(MAST_620_GFP$Gene %in% common_genes.GFP, "FDR < 0.05 in all three", "Not")

MAST_620_GFP$vplot <- case_when(MAST_620_GFP$Commonality == "FDR < 0.05 in all three" ~ "FDR < 0.05 in all three",
                                MAST_620_GFP$threshold == "FDR < 0.05" ~ "FDR < 0.05",.default = "Not Sig")


ggplot()+
  geom_point(data = MAST_620_GFP[MAST_620_GFP$vplot == "Not Sig",], aes(x = avg_log2FC, y = -log10(fdr),color=vplot))+
  geom_point(data = MAST_620_GFP[MAST_620_GFP$vplot == "FDR < 0.05",], aes(x = avg_log2FC, y = -log10(fdr),color=vplot))+
  geom_point(data = MAST_620_GFP[MAST_620_GFP$vplot == "FDR < 0.05 in all three",], aes(x = avg_log2FC, y = -log10(fdr),color=vplot))+
  scale_color_manual(values = c("red", "darkgreen", "lightgray")) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom") +
  ggtitle("A3 15dpi GFP vs dsRed Volcano Plot") +
  geom_text_repel(data = subset(MAST_620_GFP, vplot == "FDR < 0.05 in all three"),
                  aes(label = Gene,x = avg_log2FC, y = -log10(fdr)),
                  size = 3,max.overlaps = 20
  ) 



dev.off()


