#R script used to perform DGE using MAST
#
#We're using the integrated (FastMNN) .rds file 
#Additionally we also used the Wilcoxon test method, the default in Seurat, to compare with the MAST results

#Setting the seed to allow reproducibility in terms of UMAPs
set.seed(10403)

#Creating an Output.pdf file
pdf(file="DE-FastMNN-Output.pdf",compress = TRUE)

#Saving the filepath to the .RData file resulting from the main script
filePath <- c("PATH")

#Creating a directory to save the marker.csv files
OutPath <- c("PATH_II")


#Definig colors for the volcano plots
Volc.cols<-c("grey","grey","grey","red")
Volc.labels = c("NS", "NS", "NS","fdr<0.05")

#Defining variables
featuresUCOERG<- c("dsRed", "EGFP", "WPRE" )
featuresHIV<- c("gag-pol","pol","pol-vif-vpr-tat-rev","vpu-env","env","mirfp670nano","p2a-cre-ires-nef" )

#Load Packages##
library(ggplot2)
library(cli)
library(dplyr)
library(Seurat) #We used Seurat v.4.1.0 since newere version are not compatible with MAST and no bugfix has been implemented yet
library(patchwork)
library(tidyverse)
library(ggrepel)
library(gridExtra)
library(ggpmisc)
library(remotes)
library(SeuratDisk)
library(ggpubr)
library(EnhancedVolcano)
library(MAST)


#Load the .RData workspace from the previous part#
#load(filePath)

seurat.object<-readRDS(PATH)

#Create Output directory##
dir.create(OutPath)

#Differential gene expression analysis#######

##We need to add some meta.data 

##Condition (Already Implemented)
#T1/T2/T3 -> ART Treated
#A1/A2/A3 -> Acute
#U1/U2 -> Uninfected

#Condition2
#Others-> Same as Condition
#U1/U2 -> Uninfected-dsRed
#U1/U2 -> Uninfected-Unmarked



seurat.list.MNN$Condition2<-case_when(seurat.list.MNN$Count %in% c("609_count","610_count") ~ "Uninfected-Unmarked",
                                      seurat.list.MNN$Count %in% c("611_count","612_count") ~ "Uninfected-dsRed",
                                      seurat.list.MNN$MouseID %in% c("A1","A2","A3") ~ "Acute",
                                      seurat.list.MNN$MouseID %in% c("T1","T2","T3") ~ "ART Treated",)



extracted.cells <- FetchData(object = seurat.object, vars = c('orig.ident', "Count", 'nFeature_RNA', 'nCount_RNA',"Fluorescence","Treatment","SampleType","status","Condition","Condition2"), slot="counts")
extracted.cells$id <- rownames(extracted.cells)
  

#HIV+ high Acute vs. HIV- Acute
HIV.high.Acute.Cells <- extracted.cells[which(extracted.cells$Condition == "Acute" & extracted.cells$status == "HIV+ high"),]$id
seurat.object <- SetIdent(object = seurat.object, cells = HIV.high.Acute.Cells  ,value = "Acute HIV+ high")

HIV.neg.Acute.Cells <- extracted.cells[which(extracted.cells$Condition == "Acute" & extracted.cells$status == "HIV-"),]$id
seurat.object <- SetIdent(object = seurat.object, cells = HIV.neg.Acute.Cells  ,value = "Acute HIV-")



Acute.HIV.high.vs.Acute.HIV.neg.markers<-FindMarkers(seurat.object, ident.1 = "Acute HIV+ high", ident.2 = "Acute HIV-",test.use = "MAST", min.pct = 0.1, min.cells.group=1,logfc.threshold=0,  pseudocount.use= 0.001, random.seed = 12, latent.vars="nFeature_RNA")
Acute.HIV.high.vs.Acute.HIV.neg.markers = Acute.HIV.high.vs.Acute.HIV.neg.markers[ !(rownames(Acute.HIV.high.vs.Acute.HIV.neg.markers) %in% featuresHIV), ]

mousegenes<-rownames(Acute.HIV.high.vs.Acute.HIV.neg.markers)
mousegenes<-str_subset(mousegenes,pattern = "mm10*")
Acute.HIV.high.vs.Acute.HIV.neg.markers = Acute.HIV.high.vs.Acute.HIV.neg.markers[ !(rownames(Acute.HIV.high.vs.Acute.HIV.neg.markers) %in% mousegenes),]
Acute.HIV.high.vs.Acute.HIV.neg.markers$fdr <- p.adjust(Acute.HIV.high.vs.Acute.HIV.neg.markers$p_val, method = "fdr", n = length(Acute.HIV.high.vs.Acute.HIV.neg.markers$p_val))

show(EnhancedVolcano(Acute.HIV.high.vs.Acute.HIV.neg.markers, 
                     lab=rownames(Acute.HIV.high.vs.Acute.HIV.neg.markers),
                     x ="avg_log2FC", 
                     y ="p_val",subtitle = paste0("Acute HIV+ high n=",nrow(extracted.cells[which(extracted.cells$Condition == "Acute" & extracted.cells$status == "HIV+ high"),])," vs. Acute HIV- n=",nrow(extracted.cells[which(extracted.cells$Condition == "Acute" & extracted.cells$status == "HIV-"),])),
                     title = "Acute HIv+ high vs. Acute HIV-",caption = "fdr cutoff = 0.05",
                     pCutoffCol = "fdr",pCutoff=0.05,FCcutoff = 0,col = Volc.cols,cutoffLineType = "blank",legendLabels = Volc.labels))

write.csv(Acute.HIV.high.vs.Acute.HIV.neg.markers,file = paste0(OutPath,"Acute-HIV-high-vs.-Acute-HIV-neg-FastMNN-MAST-markers.csv"))

#HIV+ high Acute vs. Uninfected (dsRed + Unmarked)
HIV.high.Acute.Cells <- extracted.cells[which(extracted.cells$Condition == "Acute" & extracted.cells$status == "HIV+ high"),]$id
seurat.object <- SetIdent(object = seurat.object, cells = HIV.high.Acute.Cells  ,value = "Acute HIV+ high")

Uninfected.Cells <- extracted.cells[which(extracted.cells$Condition == "Uninfected"),]$id
seurat.object <- SetIdent(object = seurat.object, cells = Uninfected.Cells  ,value = "Uninfected")



Acute.HIV.high.vs.Uninfected.markers<-FindMarkers(seurat.object, ident.1 = "Acute HIV+ high", ident.2 = "Uninfected",test.use = "MAST", min.pct = 0.1, min.cells.group=1,logfc.threshold=0,  pseudocount.use= 0.001, random.seed = 12, latent.vars="nFeature_RNA")
Acute.HIV.high.vs.Uninfected.markers = Acute.HIV.high.vs.Uninfected.markers[ !(rownames(Acute.HIV.high.vs.Uninfected.markers) %in% featuresHIV), ]

mousegenes<-rownames(Acute.HIV.high.vs.Uninfected.markers)
mousegenes<-str_subset(mousegenes,pattern = "mm10*")
Acute.HIV.high.vs.Uninfected.markers = Acute.HIV.high.vs.Uninfected.markers[ !(rownames(Acute.HIV.high.vs.Uninfected.markers) %in% mousegenes),]
Acute.HIV.high.vs.Uninfected.markers$fdr <- p.adjust(Acute.HIV.high.vs.Uninfected.markers$p_val, method = "fdr", n = length(Acute.HIV.high.vs.Uninfected.markers$p_val))

show(EnhancedVolcano(Acute.HIV.high.vs.Uninfected.markers, 
                     lab=rownames(Acute.HIV.high.vs.Uninfected.markers),
                     x ="avg_log2FC", 
                     y ="p_val",subtitle = paste0("Acute HIV+ high n=",nrow(extracted.cells[which(extracted.cells$Condition == "Acute" & extracted.cells$status == "HIV+ high"),])," vs. Uninfected n=",nrow(extracted.cells[which(extracted.cells$Condition == "Uninfected"),])),
                     title = "Acute HIV+ high vs. Uninfected",caption = "fdr cutoff = 0.05",
                     pCutoffCol = "fdr",pCutoff=0.05,FCcutoff = 0,col = Volc.cols,cutoffLineType = "blank",legendLabels = Volc.labels))

write.csv(Acute.HIV.high.vs.Uninfected.markers,file = paste0(OutPath,"Acute-HIV-high-vs.-Uninfected-FastMNN-MAST-markers.csv"))


#Treated vs. Acute
Treated.Cells <- extracted.cells[which(extracted.cells$Condition == "ART Treated"),]$id
seurat.object <- SetIdent(object = seurat.object, cells = Treated.Cells  ,value = "ART Treated")

Acute.Cells <- extracted.cells[which(extracted.cells$Condition == "Acute"),]$id
seurat.object <- SetIdent(object = seurat.object, cells = Acute.Cells  ,value = "Acute")



Acute.vs.ARTTreated.markers<-FindMarkers(seurat.object, ident.1 = "ART Treated", ident.2 = "Acute",test.use = "MAST", min.pct = 0.1, min.cells.group=1,logfc.threshold=0,  pseudocount.use= 0.001, random.seed = 12, latent.vars="nFeature_RNA")
Acute.vs.ARTTreated.markers = Acute.vs.ARTTreated.markers[ !(rownames(Acute.vs.ARTTreated.markers) %in% featuresHIV), ]

mousegenes<-rownames(Acute.vs.ARTTreated.markers)
mousegenes<-str_subset(mousegenes,pattern = "mm10*")
Acute.vs.ARTTreated.markers = Acute.vs.ARTTreated.markers[ !(rownames(Acute.vs.ARTTreated.markers) %in% mousegenes),]
Acute.vs.ARTTreated.markers$fdr <- p.adjust(Acute.vs.ARTTreated.markers$p_val, method = "fdr", n = length(Acute.vs.ARTTreated.markers$p_val))

show(EnhancedVolcano(Acute.vs.ARTTreated.markers, 
                     lab=rownames(Acute.vs.ARTTreated.markers),
                     x ="avg_log2FC", 
                     y ="p_val",subtitle = paste0("ART Treated n=",nrow(extracted.cells[which(extracted.cells$Condition %in% c("ART Treated")),])," vs. Acute n=",nrow(extracted.cells[which(extracted.cells$Condition %in% c("Acute")),])),
                     title = "ART Treated vs. Acute",caption = "fdr cutoff = 0.05",
                     pCutoffCol = "fdr",pCutoff=0.05,FCcutoff = 0,col = Volc.cols,cutoffLineType = "blank",legendLabels = Volc.labels))

write.csv(Acute.vs.ARTTreated.markers,file = paste0(OutPath,"ART-Treated-vs.-Acute-FastMNN-MAST-markers.csv"))


#Treated vs. Uninfected dsRed

Uninfected.dsRed.Cells <- extracted.cells[which(extracted.cells$Condition2 == "Uninfected-dsRed"),]$id
seurat.object <- SetIdent(object = seurat.object, cells = Uninfected.dsRed.Cells  ,value = "Uninfected-dsRed")

Uninfected.dsRed.Cells.vs.ARTTreated.markers<-FindMarkers(seurat.object, ident.1 = "ART Treated", ident.2 = "Uninfected-dsRed",test.use = "MAST", min.pct = 0.1, min.cells.group=1,logfc.threshold=0,  pseudocount.use= 0.001, random.seed = 12, latent.vars="nFeature_RNA")
Uninfected.dsRed.Cells.vs.ARTTreated.markers = Uninfected.dsRed.Cells.vs.ARTTreated.markers[ !(rownames(Uninfected.dsRed.Cells.vs.ARTTreated.markers) %in% featuresHIV), ]

mousegenes<-rownames(Uninfected.dsRed.Cells.vs.ARTTreated.markers)
mousegenes<-str_subset(mousegenes,pattern = "mm10*")
Uninfected.dsRed.Cells.vs.ARTTreated.markers = Uninfected.dsRed.Cells.vs.ARTTreated.markers[ !(rownames(Uninfected.dsRed.Cells.vs.ARTTreated.markers) %in% mousegenes),]
Uninfected.dsRed.Cells.vs.ARTTreated.markers$fdr <- p.adjust(Uninfected.dsRed.Cells.vs.ARTTreated.markers$p_val, method = "fdr", n = length(Uninfected.dsRed.Cells.vs.ARTTreated.markers$p_val))

show(EnhancedVolcano(Uninfected.dsRed.Cells.vs.ARTTreated.markers, 
                     lab=rownames(Uninfected.dsRed.Cells.vs.ARTTreated.markers),
                     x ="avg_log2FC", 
                     y ="p_val",subtitle = paste0("ART Treated n=",nrow(extracted.cells[which(extracted.cells$Condition %in% c("ART Treated")),])," vs. Uninfected dsRed n=",nrow(extracted.cells[which(extracted.cells$Condition2 %in% c("Uninfected-dsRed")),])),
                     title = "ART Treated vs. Uninfected-dsRed",caption = "fdr cutoff = 0.05",
                     pCutoffCol = "fdr",pCutoff=0.05,FCcutoff = 0,col = Volc.cols,cutoffLineType = "blank",legendLabels = Volc.labels))

write.csv(Uninfected.dsRed.Cells.vs.ARTTreated.markers,file = paste0(OutPath,"ART-Treated-vs.-Uninfected-dsRed-FastMNN-MAST-markers.csv"))


#Treated vs. Uninfected Unmarked

Uninfected.Unmarked.Cells <- extracted.cells[which(extracted.cells$Condition2 == "Uninfected-Unmarked"),]$id
seurat.object <- SetIdent(object = seurat.object, cells = Uninfected.Unmarked.Cells  ,value = "Uninfected-Unmarked")

Uninfected.Unmarked.Cells.vs.ARTTreated.markers<-FindMarkers(seurat.object, ident.1 = "ART Treated", ident.2 = "Uninfected-Unmarked",test.use = "MAST", min.pct = 0.1, min.cells.group=1,logfc.threshold=0,  pseudocount.use= 0.001, random.seed = 12, latent.vars="nFeature_RNA")
Uninfected.Unmarked.Cells.vs.ARTTreated.markers = Uninfected.Unmarked.Cells.vs.ARTTreated.markers[ !(rownames(Uninfected.Unmarked.Cells.vs.ARTTreated.markers) %in% featuresHIV), ]

mousegenes<-rownames(Uninfected.Unmarked.Cells.vs.ARTTreated.markers)
mousegenes<-str_subset(mousegenes,pattern = "mm10*")
Uninfected.Unmarked.Cells.vs.ARTTreated.markers = Uninfected.Unmarked.Cells.vs.ARTTreated.markers[ !(rownames(Uninfected.Unmarked.Cells.vs.ARTTreated.markers) %in% mousegenes),]
Uninfected.Unmarked.Cells.vs.ARTTreated.markers$fdr <- p.adjust(Uninfected.Unmarked.Cells.vs.ARTTreated.markers$p_val, method = "fdr", n = length(Uninfected.Unmarked.Cells.vs.ARTTreated.markers$p_val))

show(EnhancedVolcano(Uninfected.Unmarked.Cells.vs.ARTTreated.markers, 
                     lab=rownames(Uninfected.Unmarked.Cells.vs.ARTTreated.markers),
                     x ="avg_log2FC", 
                     y ="p_val",subtitle = paste0("ART Treated n=",nrow(extracted.cells[which(extracted.cells$Condition %in% c("ART Treated")),])," vs. Uninfected Unmarked n=",nrow(extracted.cells[which(extracted.cells$Condition2 %in% c("Uninfected-Unmarked")),])),
                     title = "ART Treated vs. Uninfected-Unmarked",caption = "fdr cutoff = 0.05",
                     pCutoffCol = "fdr",pCutoff=0.05,FCcutoff = 0,col = Volc.cols,cutoffLineType = "blank",legendLabels = Volc.labels))

write.csv(Uninfected.Unmarked.Cells.vs.ARTTreated.markers,file = paste0(OutPath,"ART-Treated-vs.-Uninfected-Unmarked-FastMNN-MAST-markers.csv"))

#Acute vs. Uninfected dsRed

Uninfected.Unmarked.Cells.vs.Acute.markers<-FindMarkers(seurat.object, ident.1 = "Acute", ident.2 = "Uninfected-Unmarked",test.use = "MAST", min.pct = 0.1, min.cells.group=1,logfc.threshold=0,  pseudocount.use= 0.001, random.seed = 12, latent.vars="nFeature_RNA")
Uninfected.Unmarked.Cells.vs.Acute.markers = Uninfected.Unmarked.Cells.vs.Acute.markers[ !(rownames(Uninfected.Unmarked.Cells.vs.Acute.markers) %in% featuresHIV), ]

mousegenes<-rownames(Uninfected.Unmarked.Cells.vs.Acute.markers)
mousegenes<-str_subset(mousegenes,pattern = "mm10*")
Uninfected.Unmarked.Cells.vs.Acute.markers = Uninfected.Unmarked.Cells.vs.Acute.markers[ !(rownames(Uninfected.Unmarked.Cells.vs.Acute.markers) %in% mousegenes),]
Uninfected.Unmarked.Cells.vs.Acute.markers$fdr <- p.adjust(Uninfected.Unmarked.Cells.vs.Acute.markers$p_val, method = "fdr", n = length(Uninfected.Unmarked.Cells.vs.Acute.markers$p_val))

show(EnhancedVolcano(Uninfected.Unmarked.Cells.vs.Acute.markers, 
                     lab=rownames(Uninfected.Unmarked.Cells.vs.Acute.markers),
                     x ="avg_log2FC", 
                     y ="p_val",subtitle = paste0("Acute n=",nrow(extracted.cells[which(extracted.cells$Condition %in% c("Acute")),])," vs. Uninfected Unmarked n=",nrow(extracted.cells[which(extracted.cells$Condition2 %in% c("Uninfected-Unmarked")),])),
                     title = "Acute vs. Uninfected-Unmarked",caption = "fdr cutoff = 0.05",
                     pCutoffCol = "fdr",pCutoff=0.05,FCcutoff = 0,col = Volc.cols,cutoffLineType = "blank",legendLabels = Volc.labels))

write.csv(Uninfected.Unmarked.Cells.vs.Acute.markers,file = paste0(OutPath,"Acute-vs.-Uninfected-Unmarked-FastMNN-MAST-markers.csv"))


#Acute vs. Uninfected Unmarked

Uninfected.dsRed.Cells.vs.Acute.markers<-FindMarkers(seurat.object, ident.1 = "Acute", ident.2 = "Uninfected-dsRed",test.use = "MAST", min.pct = 0.1, min.cells.group=1,logfc.threshold=0,  pseudocount.use= 0.001, random.seed = 12, latent.vars="nFeature_RNA")
Uninfected.dsRed.Cells.vs.Acute.markers = Uninfected.dsRed.Cells.vs.Acute.markers[ !(rownames(Uninfected.dsRed.Cells.vs.Acute.markers) %in% featuresHIV), ]

mousegenes<-rownames(Uninfected.dsRed.Cells.vs.Acute.markers)
mousegenes<-str_subset(mousegenes,pattern = "mm10*")
Uninfected.dsRed.Cells.vs.Acute.markers = Uninfected.dsRed.Cells.vs.Acute.markers[ !(rownames(Uninfected.dsRed.Cells.vs.Acute.markers) %in% mousegenes),]
Uninfected.dsRed.Cells.vs.Acute.markers$fdr <- p.adjust(Uninfected.dsRed.Cells.vs.Acute.markers$p_val, method = "fdr", n = length(Uninfected.dsRed.Cells.vs.Acute.markers$p_val))

show(EnhancedVolcano(Uninfected.dsRed.Cells.vs.Acute.markers, 
                     lab=rownames(Uninfected.dsRed.Cells.vs.Acute.markers),
                     x ="avg_log2FC", 
                     y ="p_val",subtitle = paste0("Acute n=",nrow(extracted.cells[which(extracted.cells$Condition %in% c("Acute")),])," vs. Uninfected dsRed n=",nrow(extracted.cells[which(extracted.cells$Condition2 %in% c("Uninfected-dsRed")),])),
                     title = "Acute vs. Uninfected-dsRed",caption = "fdr cutoff = 0.05",
                     pCutoffCol = "fdr",pCutoff=0.05,FCcutoff = 0,col = Volc.cols,cutoffLineType = "blank",legendLabels = Volc.labels))

write.csv(Uninfected.dsRed.Cells.vs.Acute.markers,file = paste0(OutPath,"Acute-vs.-Uninfected-dsRed-FastMNN-MAST-markers.csv"))

#HIV+ high vs Uninfected Unmarked
HIVhigh.Acute.Cells <- extracted.cells[which(extracted.cells$Condition2 == "Acute" & extracted.cells$status == "HIV+ high" ),]$id
seurat.object <- SetIdent(object = seurat.object, cells = HIVhigh.Acute.Cells  ,value = "HIV+ high Acute")

Uninfected.Unmarked.Cells.vs.HIVhighAcute.markers<-FindMarkers(seurat.object, ident.1 = "HIV+ high Acute", ident.2 = "Uninfected-Unmarked",test.use = "MAST", min.pct = 0.1, min.cells.group=1,logfc.threshold=0,  pseudocount.use= 0.001, random.seed = 12, latent.vars="nFeature_RNA")
Uninfected.Unmarked.Cells.vs.HIVhighAcute.markers = Uninfected.Unmarked.Cells.vs.HIVhighAcute.markers[ !(rownames(Uninfected.Unmarked.Cells.vs.HIVhighAcute.markers) %in% featuresHIV), ]

mousegenes<-rownames(Uninfected.Unmarked.Cells.vs.HIVhighAcute.markers)
mousegenes<-str_subset(mousegenes,pattern = "mm10*")
Uninfected.Unmarked.Cells.vs.HIVhighAcute.markers = Uninfected.Unmarked.Cells.vs.HIVhighAcute.markers[ !(rownames(Uninfected.Unmarked.Cells.vs.HIVhighAcute.markers) %in% mousegenes),]
Uninfected.Unmarked.Cells.vs.HIVhighAcute.markers$fdr <- p.adjust(Uninfected.Unmarked.Cells.vs.HIVhighAcute.markers$p_val, method = "fdr", n = length(Uninfected.Unmarked.Cells.vs.HIVhighAcute.markers$p_val))

show(EnhancedVolcano(Uninfected.Unmarked.Cells.vs.HIVhighAcute.markers, 
                     lab=rownames(Uninfected.Unmarked.Cells.vs.HIVhighAcute.markers),
                     x ="avg_log2FC", 
                     y ="p_val",subtitle = paste0("HIV+ high Acute n=",nrow(extracted.cells[which(extracted.cells$Condition2 == "Acute" & extracted.cells$status == "HIV+ high"),])," vs. Uninfected Unmarked n=",nrow(extracted.cells[which(extracted.cells$Condition2 %in% c("Uninfected-Unmarked")),])),
                     title = "HIV+ high Acute vs. Uninfected-Unmarked",caption = "fdr cutoff = 0.05",
                     pCutoffCol = "fdr",pCutoff=0.05,FCcutoff = 0,col = Volc.cols,cutoffLineType = "blank",legendLabels = Volc.labels))

write.csv(Uninfected.Unmarked.Cells.vs.HIVhighAcute.markers,file = paste0(OutPath,"/HIVhighAcute-vs.-Uninfected-Unmarked-FastMNN-MAST-markers.csv"))

#HIV+ high vs Uninfected dsRed
HIVhigh.Acute.Cells <- extracted.cells[which(extracted.cells$Condition2 == "Acute" & extracted.cells$status == "HIV+ high" ),]$id
seurat.object <- SetIdent(object = seurat.object, cells = HIVhigh.Acute.Cells  ,value = "HIV+ high Acute")

Uninfected.dsRed.Cells.vs.HIVhighAcute.markers<-FindMarkers(seurat.object, ident.1 = "HIV+ high Acute", ident.2 = "Uninfected-dsRed",test.use = "MAST", min.pct = 0.1, min.cells.group=1,logfc.threshold=0,  pseudocount.use= 0.001, random.seed = 12, latent.vars="nFeature_RNA")
Uninfected.dsRed.Cells.vs.HIVhighAcute.markers = Uninfected.dsRed.Cells.vs.HIVhighAcute.markers[ !(rownames(Uninfected.dsRed.Cells.vs.HIVhighAcute.markers) %in% featuresHIV), ]

mousegenes<-rownames(Uninfected.dsRed.Cells.vs.HIVhighAcute.markers)
mousegenes<-str_subset(mousegenes,pattern = "mm10*")
Uninfected.dsRed.Cells.vs.HIVhighAcute.markers = Uninfected.dsRed.Cells.vs.HIVhighAcute.markers[ !(rownames(Uninfected.dsRed.Cells.vs.HIVhighAcute.markers) %in% mousegenes),]
Uninfected.dsRed.Cells.vs.HIVhighAcute.markers$fdr <- p.adjust(Uninfected.dsRed.Cells.vs.HIVhighAcute.markers$p_val, method = "fdr", n = length(Uninfected.dsRed.Cells.vs.HIVhighAcute.markers$p_val))

show(EnhancedVolcano(Uninfected.dsRed.Cells.vs.HIVhighAcute.markers, 
                     lab=rownames(Uninfected.dsRed.Cells.vs.HIVhighAcute.markers),
                     x ="avg_log2FC", 
                     y ="p_val",subtitle = paste0("HIV+ high Acute n=",nrow(extracted.cells[which(extracted.cells$Condition2 == "Acute" & extracted.cells$status == "HIV+ high"),])," vs. Uninfected dsRed n=",nrow(extracted.cells[which(extracted.cells$Condition2 %in% c("Uninfected-dsRed")),])),
                     title = "HIV+ high Acute vs. Uninfected-dsRed",caption = "fdr cutoff = 0.05",
                     pCutoffCol = "fdr",pCutoff=0.05,FCcutoff = 0,col = Volc.cols,cutoffLineType = "blank",legendLabels = Volc.labels))

write.csv(Uninfected.dsRed.Cells.vs.HIVhighAcute.markers,file = paste0(OutPath,"/HIVhighAcute-vs.-Uninfected-dsRed-FastMNN-MAST-markers.csv"))


#Save of the R workspace#####
save.image("MAST-DE-FastMNN.RData")

#Finish pdf###
dev.off()

#Acute vs. Uninfected Unmarked
#We make a copy with only those that are fdr<0.05
Acute.vs..Uninfected.Unmarked.FastMNN.MAST.markers$threshold <- ifelse(Acute.vs..Uninfected.Unmarked.FastMNN.MAST.markers$fdr < 0.05, "FDR < 0.05", "Not Sig")

#
Acute.vs..Uninfected.Unmarked.FastMNN.MAST.markers_sig<-subset.data.frame(Acute.vs..Uninfected.Unmarked.FastMNN.MAST.markers, threshold == "FDR < 0.05")

#We make copies and only take the genes columns
Acute.vs..Uninfected.Unmarked.FastMNN.MAST.markers_G<-subset.data.frame(Acute.vs..Uninfected.Unmarked.FastMNN.MAST.markers_sig,select = "X")


#VolcanoPlot
ggplot()+
  geom_point(data = Acute.vs..Uninfected.Unmarked.FastMNN.MAST.markers, aes(x = avg_log2FC, y = -log10(fdr),color=threshold))+
  scale_color_manual(values = c("red", "lightgray")) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom") +
  ggtitle("Acute vs Uninfected-Unmarked FastMNN Volcano Plot") +
  geom_text_repel(data = subset(Acute.vs..Uninfected.Unmarked.FastMNN.MAST.markers, threshold == "FDR < 0.05"),
                  aes(label = X,x = avg_log2FC, y = -log10(fdr)),
                  size = 3,max.overlaps = 20) 


#Acute vs. Uninfected dsRed
#We make a copy with only those that are fdr<0.05
Acute.vs..Uninfected.dsRed.FastMNN.MAST.markers$threshold <- ifelse(Acute.vs..Uninfected.dsRed.FastMNN.MAST.markers$fdr < 0.05, "FDR < 0.05", "Not Sig")

#
Acute.vs..Uninfected.dsRed.FastMNN.MAST.markers_sig<-subset.data.frame(Acute.vs..Uninfected.dsRed.FastMNN.MAST.markers, threshold == "FDR < 0.05")

#We make copies and only take the genes columns
Acute.vs..Uninfected.dsRed.FastMNN.MAST.markers_G<-subset.data.frame(Acute.vs..Uninfected.dsRed.FastMNN.MAST.markers_sig,select = "X")


#VolcanoPlot
ggplot()+
  geom_point(data = Acute.vs..Uninfected.dsRed.FastMNN.MAST.markers, aes(x = avg_log2FC, y = -log10(fdr),color=threshold))+
  scale_color_manual(values = c("red", "lightgray")) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom") +
  ggtitle("Acute vs Uninfected-dsRed FastMNN Volcano Plot") +
  geom_text_repel(data = subset(Acute.vs..Uninfected.dsRed.FastMNN.MAST.markers, threshold == "FDR < 0.05"),
                  aes(label = X,x = avg_log2FC, y = -log10(fdr)),
                  size = 3,max.overlaps = 20) 



#ART-Treated-vs.-Acute-FastMNN-MAST
#We make a copy with only those that are fdr<0.05
ART.Treated.vs..Acute.FastMNN.MAST.markers$threshold <- ifelse(ART.Treated.vs..Acute.FastMNN.MAST.markers$fdr < 0.05, "FDR < 0.05", "Not Sig")

#
ART.Treated.vs..Acute.FastMNN.MAST.markers_sig<-subset.data.frame(ART.Treated.vs..Acute.FastMNN.MAST.markers, threshold == "FDR < 0.05")

#We make copies and only take the genes columns
ART.Treated.vs..Acute.FastMNN.MAST.markers_G<-subset.data.frame(ART.Treated.vs..Acute.FastMNN.MAST.markers_sig,select = "X")


#VolcanoPlot
ggplot()+
  geom_point(data = ART.Treated.vs..Acute.FastMNN.MAST.markers, aes(x = avg_log2FC, y = -log10(fdr),color=threshold))+
  scale_color_manual(values = c("red", "lightgray")) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom") +
  ggtitle("ART.Treated.vs..Acute.FastMNN FastMNN Volcano Plot") +
  geom_text_repel(data = subset(ART.Treated.vs..Acute.FastMNN.MAST.markers, threshold == "FDR < 0.05"),
                  aes(label = X,x = avg_log2FC, y = -log10(fdr)),
                  size = 3,max.overlaps = 20) 


#ART-Treated-vs.-Uninfected dsRed-FastMNN-MAST
#We make a copy with only those that are fdr<0.05
ART.Treated.vs..Uninfected.dsRed.FastMNN.MAST.markers$threshold <- ifelse(ART.Treated.vs..Uninfected.dsRed.FastMNN.MAST.markers$fdr < 0.05, "FDR < 0.05", "Not Sig")

#
ART.Treated.vs..Uninfected.dsRed.FastMNN.MAST.markers_sig<-subset.data.frame(ART.Treated.vs..Uninfected.dsRed.FastMNN.MAST.markers, threshold == "FDR < 0.05")

#We make copies and only take the genes columns
ART.Treated.vs..Uninfected.dsRed.FastMNN.MAST.markers_G<-subset.data.frame(ART.Treated.vs..Uninfected.dsRed.FastMNN.MAST.markers_sig,select = "X")


#VolcanoPlot
ggplot()+
  geom_point(data = ART.Treated.vs..Uninfected.dsRed.FastMNN.MAST.markers, aes(x = avg_log2FC, y = -log10(fdr),color=threshold))+
  scale_color_manual(values = c("red", "lightgray")) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom") +
  ggtitle("ART.Treated.vs..Uninfected.dsRed.FastMNN FastMNN Volcano Plot") +
  geom_text_repel(data = subset(ART.Treated.vs..Uninfected.dsRed.FastMNN.MAST.markers, threshold == "FDR < 0.05"),
                  aes(label = X,x = avg_log2FC, y = -log10(fdr)),
                  size = 3,max.overlaps = 20) 



#ART-Treated-vs.-Uninfected Unmarked-FastMNN-MAST
#We make a copy with only those that are fdr<0.05
ART.Treated.vs..Uninfected.Unmarked.FastMNN.MAST.markers$threshold <- ifelse(ART.Treated.vs..Uninfected.Unmarked.FastMNN.MAST.markers$fdr < 0.05, "FDR < 0.05", "Not Sig")

#
ART.Treated.vs..Uninfected.Unmarked.FastMNN.MAST.markers_sig<-subset.data.frame(ART.Treated.vs..Uninfected.Unmarked.FastMNN.MAST.markers, threshold == "FDR < 0.05")

#We make copies and only take the genes columns
ART.Treated.vs..Uninfected.Unmarked.FastMNN.MAST.markers_G<-subset.data.frame(ART.Treated.vs..Uninfected.Unmarked.FastMNN.MAST.markers_sig,select = "X")


#VolcanoPlot
ggplot()+
  geom_point(data = ART.Treated.vs..Uninfected.Unmarked.FastMNN.MAST.markers, aes(x = avg_log2FC, y = -log10(fdr),color=threshold))+
  scale_color_manual(values = c("red", "lightgray")) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom") +
  ggtitle("ART.Treated.vs..Uninfected.Unmarked.FastMNN FastMNN Volcano Plot") +
  geom_text_repel(data = subset(ART.Treated.vs..Uninfected.Unmarked.FastMNN.MAST.markers, threshold == "FDR < 0.05"),
                  aes(label = X,x = avg_log2FC, y = -log10(fdr)),
                  size = 3,max.overlaps = 20) 


##END OF SCRIPT####



