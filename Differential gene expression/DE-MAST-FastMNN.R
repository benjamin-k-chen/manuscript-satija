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



extracted.cells <- FetchData(object = seurat.list.MNN, vars = c("Count",'MouseID',"Fluorescence","Treatment","SampleType","status","Condition","Condition2",'predicted.celltype.l2','seurat_clusters'), slot="counts")
extracted.cells$id <- rownames(extracted.cells)
  

#HIV+ high vs. HIV- (within Acute)

#Extract cells of each predicted.celltype.l2 from an Acute Mouse (HIV+high (a) and HIV- (b))
  a<-extracted.cells[which(extracted.cells$Condition == "Acute" & extracted.cells$status == "HIV+ high"),]$id
  b<-extracted.cells[which(extracted.cells$Condition == "Acute" & extracted.cells$status == "HIV-"),]$id
  
  assign(paste0(unique(extracted.cells$Condition),"-HIV.high.Acute.Cells"),a)
  assign(paste0(unique(extracted.cells$Condition),"-HIV.neg.Acute.Cells"),b)
  
  print('Acute HIV+ high vs. HIV- - performing DEG ')
    
  #We perform the MAST DEG analysis on a copy of the seurat.object
  #We set cells a and b as identities to compare
  seurat.object <- SetIdent(object = seurat.list.MNN, cells = a  ,value = "HIV+ high")
  seurat.object <- SetIdent(object = seurat.object, cells = b  ,value = "HIV-")
    
    
  markers<-FindMarkers(seurat.object, ident.1 = "HIV+ high", ident.2 = "HIV-",test.use = "MAST", min.pct = 0.1, min.cells.group=1,logfc.threshold=0,  pseudocount.use= 0.001, random.seed = 12, latent.vars="nFeature_RNA")
  markers = markers[ !(rownames(markers) %in% featuresHIV),]
  
  mousegenes<-rownames(markers)
  mousegenes<-str_subset(mousegenes,pattern = "mm10*")
  markers = markers[ !(rownames(markers) %in% mousegenes),]
  markers$fdr <- p.adjust(markers$p_val, method = "fdr", n = length(markers$p_val))
    
  show(EnhancedVolcano(markers, 
                       lab=rownames(markers),
                       x ="avg_log2FC", 
                       y ="p_val",subtitle = paste0("HIV+ high n=",length(a)," vs. HIV- n=",length(b)),
                       title = 'Acute Mice',caption = "fdr cutoff = 0.05",
                       pCutoffCol = "fdr",pCutoff=0.05,FCcutoff = 0,col = Volc.cols,cutoffLineType = "blank",legendLabels = Volc.labels))
  
  write.csv(markers,file = paste0(OutPath,'/',"All_Acute_HIVhighpos_vs_HIVneg_Acute.csv"))

  #HIV+ high vs. HIV- (Uninf. dsRed)
  #Extract cells of each predicted.celltype.l2 from an Acute Mouse (HIV+high (a) and HIV- (b))
  a<-extracted.cells[which(extracted.cells$Condition == "Acute" & extracted.cells$status == "HIV+ high"),]$id
  b<-extracted.cells[which(extracted.cells$Condition2 == "Uninfected-dsRed" & extracted.cells$status == "HIV-"),]$id
  
  assign(paste0(unique(extracted.cells$Condition),"-HIV.high.Acute.Cells"),a)
  assign(paste0(unique(extracted.cells$Condition2),"-HIV.neg.Cells"),b)
  
  print('Acute HIV+ high vs. HIV- Uninfected-dsRed - performing DEG ')
  
  #We perform the MAST DEG analysis on a copy of the seurat.object
  #We set cells a and b as identities to compare
  seurat.object <- SetIdent(object = seurat.list.MNN, cells = a  ,value = "HIV+ high")
  seurat.object <- SetIdent(object = seurat.object, cells = b  ,value = "HIV-")
  
  
  markers<-FindMarkers(seurat.object, ident.1 = "HIV+ high", ident.2 = "HIV-",test.use = "MAST", min.pct = 0.1, min.cells.group=1,logfc.threshold=0,  pseudocount.use= 0.001, random.seed = 12, latent.vars="nFeature_RNA")
  markers = markers[ !(rownames(markers) %in% featuresHIV),]
  
  mousegenes<-rownames(markers)
  mousegenes<-str_subset(mousegenes,pattern = "mm10*")
  markers = markers[ !(rownames(markers) %in% mousegenes),] 
  markers$fdr <- p.adjust(markers$p_val, method = "fdr", n = length(markers$p_val))
  
  show(EnhancedVolcano(markers, 
                       lab=rownames(markers),
                       x ="avg_log2FC", 
                       y ="p_val",subtitle = paste0("HIV+ high n=",length(a)," vs. HIV- n=",length(b)),
                       title = 'Acute HIV+ high vs. Uninf. dsRed',caption = "fdr cutoff = 0.05",
                       pCutoffCol = "fdr",pCutoff=0.05,FCcutoff = 0,col = Volc.cols,cutoffLineType = "blank",legendLabels = Volc.labels))
  
  write.csv(markers,file = paste0(OutPath,'/',"All_Acute_HIVhighpos_vs_HIVneg_Uninf_dsRed.csv"))


#HIV+ high vs. HIV- (within Acute) for each predicted.celltype.l2
for(i in 1:length(extracted.cells.split.celltype)){
  #Extract cells of each predicted.celltype.l2 from an Acute Mouse (HIV+high (a) and HIV- (b))
  a<-extracted.cells.split.celltype[[i]][which(extracted.cells.split.celltype[[i]]$Condition == "Acute" & extracted.cells.split.celltype[[i]]$status == "HIV+ high"),]$id
  b<-extracted.cells.split.celltype[[i]][which(extracted.cells.split.celltype[[i]]$Condition == "Acute" & extracted.cells.split.celltype[[i]]$status == "HIV-"),]$id
  
  assign(paste0(unique(extracted.cells.split.celltype[[i]]$predicted.celltype.l2),"-HIV.high.Acute.Cells"),a)
  print(a)
  
  if(length(a)>1){
    print(paste0(unique(extracted.cells.split.celltype[[i]]$predicted.celltype.l2),' - performing DEG with Acute'))
    
    #We perform the MAST DEG analysis on a copy of the seurat.object
    #We set cells a and b as identities to compare
    seurat.object <- SetIdent(object = seurat.list.MNN, cells = a  ,value = "HIV+ high")
    
    seurat.object <- SetIdent(object = seurat.object, cells = b  ,value = "HIV-")
    
  
    markers<-FindMarkers(seurat.object, ident.1 = "HIV+ high", ident.2 = "HIV-",test.use = "MAST", min.pct = 0.1, min.cells.group=1,logfc.threshold=0,  pseudocount.use= 0.001, random.seed = 12, latent.vars="nFeature_RNA")
    markers = markers[ !(rownames(markers) %in% featuresHIV),]
    
    mousegenes<-rownames(markers)
    mousegenes<-str_subset(mousegenes,pattern = "mm10*")
    markers = markers[ !(rownames(markers) %in% mousegenes),]
    markers$fdr <- p.adjust(markers$p_val, method = "fdr", n = length(markers$p_val))
    
    show(EnhancedVolcano(markers, 
                         lab=rownames(markers),
                         x ="avg_log2FC", 
                         y ="p_val",subtitle = paste0("Acute HIV+ high n=",length(a)," vs. Acute HIV- n=",length(b)),
                         title = paste0(unique(extracted.cells.split.celltype[[i]]$predicted.celltype.l2)),caption = "fdr cutoff = 0.05",
                         pCutoffCol = "fdr",pCutoff=0.05,FCcutoff = 0,col = Volc.cols,cutoffLineType = "blank",legendLabels = Volc.labels))
    
    write.csv(markers,file = paste0(OutPath,'/',unique(extracted.cells.split.celltype[[i]]$predicted.celltype.l2),"_HIVhighpos_vs_HIVneg_Acute.csv"))
  }
}

#HIV+ high vs. HIV- (regardless of condition) for each predicted.celltype.l2
for(i in 1:length(extracted.cells.split.celltype)){
    #Extract cells of each predicted.celltype.l2 from an Acute Mouse (HIV+high (a) and HIV- (b))
    a<-extracted.cells.split.celltype[[i]][which(extracted.cells.split.celltype[[i]]$Condition == "Acute" & extracted.cells.split.celltype[[i]]$status == "HIV+ high"),]$id
    b<-extracted.cells.split.celltype[[i]][which(extracted.cells.split.celltype[[i]]$status == "HIV-"),]$id
    
    assign(paste0(unique(extracted.cells.split.celltype[[i]]$predicted.celltype.l2),"-HIV.high.Acute.Cells"),a)
    
    if(length(a)>1){
      print(paste0(unique(extracted.cells.split.celltype[[i]]$predicted.celltype.l2),'- performing DEG all cells in cell type'))
      
      #We perform the MAST DEG analysis on a copy of the seurat.object
      #We set cells a and b as identities to compare
      seurat.object <- SetIdent(object = seurat.list.MNN, cells = a  ,value = "HIV+ high")
      
      seurat.object <- SetIdent(object = seurat.object, cells = b  ,value = "HIV-")
      
      
      markers<-FindMarkers(seurat.object, ident.1 = "HIV+ high", ident.2 = "HIV-",test.use = "MAST", min.pct = 0.1, min.cells.group=1,logfc.threshold=0,  pseudocount.use= 0.001, random.seed = 12, latent.vars="nFeature_RNA")
      markers = markers[ !(rownames(markers) %in% featuresHIV),]
      
      mousegenes<-rownames(markers)
      mousegenes<-str_subset(mousegenes,pattern = "mm10*")
      markers = markers[ !(rownames(markers) %in% mousegenes),]
      markers$fdr <- p.adjust(markers$p_val, method = "fdr", n = length(markers$p_val))
      
      show(EnhancedVolcano(markers, 
                           lab=rownames(markers),
                           x ="avg_log2FC", 
                           y ="p_val",subtitle = paste0("Acute HIV+ high n=",length(a)," vs. HIV- n=",length(b)),
                           title = paste0(unique(extracted.cells.split.celltype[[i]]$predicted.celltype.l2)),caption = "fdr cutoff = 0.05",
                           pCutoffCol = "fdr",pCutoff=0.05,FCcutoff = 0,col = Volc.cols,cutoffLineType = "blank",legendLabels = Volc.labels))
      
      write.csv(markers,file = paste0(OutPath,'/',unique(extracted.cells.split.celltype[[i]]$predicted.celltype.l2),"_HIVhighpos_vs_HIVneg_all.csv"))
    }
  }

#HIV+ high vs. HIV- (Uninf. dsRed)
for(i in 1:length(extracted.cells.split.celltype)){
    #Extract cells of each predicted.celltype.l2 from an Acute Mouse (HIV+high (a) and HIV- (b))
    a<-extracted.cells.split.celltype[[i]][which(extracted.cells.split.celltype[[i]]$Condition == "Acute" & extracted.cells.split.celltype[[i]]$status == "HIV+ high"),]$id
    b<-extracted.cells.split.celltype[[i]][which(extracted.cells.split.celltype[[i]]$Condition2 == "Uninfected-dsRed" & extracted.cells.split.celltype[[i]]$status == "HIV-"),]$id
    
    assign(paste0(unique(extracted.cells.split.celltype[[i]]$predicted.celltype.l2),"-HIV.high.Acute.Cells"),a)
    
    if(length(a)>1){
      print(paste0(unique(extracted.cells.split.celltype[[i]]$predicted.celltype.l2),' - performing DEG with Uninfected dsRed'))
      
      #We perform the MAST DEG analysis on a copy of the seurat.object
      #We set cells a and b as identities to compare
      seurat.object <- SetIdent(object = seurat.list.MNN, cells = a  ,value = "HIV+ high")
      
      seurat.object <- SetIdent(object = seurat.object, cells = b  ,value = "HIV-")
      
      
      markers<-FindMarkers(seurat.object, ident.1 = "HIV+ high", ident.2 = "HIV-",test.use = "MAST", min.pct = 0.1, min.cells.group=1,logfc.threshold=0,  pseudocount.use= 0.001, random.seed = 12, latent.vars="nFeature_RNA")
      markers = markers[ !(rownames(markers) %in% featuresHIV),]
      
      mousegenes<-rownames(markers)
      mousegenes<-str_subset(mousegenes,pattern = "mm10*")
      markers = markers[ !(rownames(markers) %in% mousegenes),]
      markers$fdr <- p.adjust(markers$p_val, method = "fdr", n = length(markers$p_val))
      
      show(EnhancedVolcano(markers, 
                           lab=rownames(markers),
                           x ="avg_log2FC", 
                           y ="p_val",subtitle = paste0("Acute HIV+ high n=",length(a)," vs. Uninf. dsRed HIV- n=",length(b)),
                           title = paste0(unique(extracted.cells.split.celltype[[i]]$predicted.celltype.l2)),caption = "fdr cutoff = 0.05",
                           pCutoffCol = "fdr",pCutoff=0.05,FCcutoff = 0,col = Volc.cols,cutoffLineType = "blank",legendLabels = Volc.labels))
      
      write.csv(markers,file = paste0(OutPath,'/',unique(extracted.cells.split.celltype[[i]]$predicted.celltype.l2),"_HIVhighpos_vs_HIVneg_Acute_Uninf_dsRed.csv"))
    }
  }

######  

#HIV+ high vs. HIV- (within Acute) for each seurat_cluster
#First we need to resplit extracted.cells by seurat.cluster
extracted.cells.split.cluster<-split(extracted.cells,extracted.cells$seurat_clusters)

#HIV+ high vs. HIV- (within Acute) for each seurat_cluster
for(i in 1:length(extracted.cells.split.cluster)){
  #Extract cells of each seurat_clusters from an Acute Mouse (HIV+high (a) and HIV- (b))
  a<-extracted.cells.split.cluster[[i]][which(extracted.cells.split.cluster[[i]]$Condition == "Acute" & extracted.cells.split.cluster[[i]]$status == "HIV+ high"),]$id
  b<-extracted.cells.split.cluster[[i]][which(extracted.cells.split.cluster[[i]]$Condition == "Acute" & extracted.cells.split.cluster[[i]]$status == "HIV-"),]$id
  
  assign(paste0('Cluster_',unique(extracted.cells.split.cluster[[i]]$seurat_clusters),"-HIV.high.Acute.Cells"),a)
  
  if(length(a)>1){
    print(paste0('Cluster',unique(extracted.cells.split.cluster[[i]]$seurat_clusters),' - performing DEG with Acute'))
    
    #We perform the MAST DEG analysis on a copy of the seurat.object
    #We set cells a and b as identities to compare
    seurat.object <- SetIdent(object = seurat.list.MNN, cells = a  ,value = "HIV+ high")
    
    seurat.object <- SetIdent(object = seurat.object, cells = b  ,value = "HIV-")
    
    
    markers<-FindMarkers(seurat.object, ident.1 = "HIV+ high", ident.2 = "HIV-",test.use = "MAST", min.pct = 0.1, min.cells.group=1,logfc.threshold=0,  pseudocount.use= 0.001, random.seed = 12, latent.vars="nFeature_RNA")
    markers = markers[ !(rownames(markers) %in% featuresHIV),]
    
    mousegenes<-rownames(markers)
    mousegenes<-str_subset(mousegenes,pattern = "mm10*")
    markers = markers[ !(rownames(markers) %in% mousegenes),]
    markers$fdr <- p.adjust(markers$p_val, method = "fdr", n = length(markers$p_val))
    
    show(EnhancedVolcano(markers, 
                         lab=rownames(markers),
                         x ="avg_log2FC", 
                         y ="p_val",subtitle = paste0("Acute HIV+ high n=",length(a)," vs. Acute HIV- n=",length(b)),
                         title = paste0('Cluster ',unique(extracted.cells.split.cluster[[i]]$seurat_clusters)),caption = "fdr cutoff = 0.05",
                         pCutoffCol = "fdr",pCutoff=0.05,FCcutoff = 0,col = Volc.cols,cutoffLineType = "blank",legendLabels = Volc.labels))
    
    write.csv(markers,file = paste0(OutPath,'/','Cluster_',unique(extracted.cells.split.cluster[[i]]$seurat_clusters),"_HIVhighpos_vs_HIVneg_Acute.csv"))
  }
}

#HIV+ high vs. HIV- (regardless of condition) for each seurat_cluster
for(i in 1:length(extracted.cells.split.cluster)){
    #Extract cells of each seurat_clusters from an Acute Mouse (HIV+high (a) and HIV- (b))
    a<-extracted.cells.split.cluster[[i]][which(extracted.cells.split.cluster[[i]]$Condition == "Acute" & extracted.cells.split.cluster[[i]]$status == "HIV+ high"),]$id
    b<-extracted.cells.split.cluster[[i]][which(extracted.cells.split.cluster[[i]]$status == "HIV-"),]$id
    
    assign(paste0('Cluster_',unique(extracted.cells.split.cluster[[i]]$seurat_clusters),"-HIV.high.Acute.Cells"),a)
    
    if(length(a)>1){
      print(paste0('Cluster',unique(extracted.cells.split.cluster[[i]]$seurat_clusters),' - performing DEG all cells in cluster'))
      
      #We perform the MAST DEG analysis on a copy of the seurat.object
      #We set cells a and b as identities to compare
      seurat.object <- SetIdent(object = seurat.list.MNN, cells = a  ,value = "HIV+ high")
      
      seurat.object <- SetIdent(object = seurat.object, cells = b  ,value = "HIV-")
      
      
      markers<-FindMarkers(seurat.object, ident.1 = "HIV+ high", ident.2 = "HIV-",test.use = "MAST", min.pct = 0.1, min.cells.group=1,logfc.threshold=0,  pseudocount.use= 0.001, random.seed = 12, latent.vars="nFeature_RNA")
      markers = markers[ !(rownames(markers) %in% featuresHIV),]
      
      mousegenes<-rownames(markers)
      mousegenes<-str_subset(mousegenes,pattern = "mm10*")
      markers = markers[ !(rownames(markers) %in% mousegenes),]
      markers$fdr <- p.adjust(markers$p_val, method = "fdr", n = length(markers$p_val))
      
      show(EnhancedVolcano(markers, 
                           lab=rownames(markers),
                           x ="avg_log2FC", 
                           y ="p_val",subtitle = paste0("Acute HIV+ high n=",length(a)," vs. HIV- n=",length(b)),
                           title = paste0('Cluster ',unique(extracted.cells.split.cluster[[i]]$seurat_clusters)),caption = "fdr cutoff = 0.05",
                           pCutoffCol = "fdr",pCutoff=0.05,FCcutoff = 0,col = Volc.cols,cutoffLineType = "blank",legendLabels = Volc.labels))
      
      write.csv(markers,file = paste0(OutPath,'/','Cluster_',unique(extracted.cells.split.cluster[[i]]$seurat_clusters),"_HIVhighpos_vs_HIVneg_all.csv"))
    }
  }
  
#HIV+ high vs. HIV- (Uninf. dsRed) for each seurat_cluster
for(i in 1:length(extracted.cells.split.cluster)){
  #Extract cells of each seurat_clusters from an Acute Mouse (HIV+high (a) and HIV- (b))
  a<-extracted.cells.split.cluster[[i]][which(extracted.cells.split.cluster[[i]]$Condition == "Acute" & extracted.cells.split.cluster[[i]]$status == "HIV+ high"),]$id
  b<-extracted.cells.split.cluster[[i]][which(extracted.cells.split.cluster[[i]]$Condition2 == "Uninfected-dsRed" & extracted.cells.split.cluster[[i]]$status == "HIV-"),]$id
  
  assign(paste0('Cluster_',unique(extracted.cells.split.cluster[[i]]$seurat_clusters),"-HIV.high.Acute.Cells"),a)
  
  if(length(a)>1){
    print(paste0('Cluster',unique(extracted.cells.split.cluster[[i]]$seurat_clusters),' - performing DEG with Uninfected dsRed'))
    
    #We perform the MAST DEG analysis on a copy of the seurat.object
    #We set cells a and b as identities to compare
    seurat.object <- SetIdent(object = seurat.list.MNN, cells = a  ,value = "HIV+ high")
    
    seurat.object <- SetIdent(object = seurat.object, cells = b  ,value = "HIV-")
    
    
    markers<-FindMarkers(seurat.object, ident.1 = "HIV+ high", ident.2 = "HIV-",test.use = "MAST", min.pct = 0.1, min.cells.group=1,logfc.threshold=0,  pseudocount.use= 0.001, random.seed = 12, latent.vars="nFeature_RNA")
    markers = markers[ !(rownames(markers) %in% featuresHIV),]
    
    mousegenes<-rownames(markers)
    mousegenes<-str_subset(mousegenes,pattern = "mm10*")
    markers = markers[ !(rownames(markers) %in% mousegenes),]
    markers$fdr <- p.adjust(markers$p_val, method = "fdr", n = length(markers$p_val))
    
    show(EnhancedVolcano(markers, 
                         lab=rownames(markers),
                         x ="avg_log2FC", 
                         y ="p_val",subtitle = paste0("Acute HIV+ high n=",length(a)," vs. Uninfected dsRed HIV- n=",length(b)),
                         title = paste0('Cluster ',unique(extracted.cells.split.cluster[[i]]$seurat_clusters)),caption = "fdr cutoff = 0.05",
                         pCutoffCol = "fdr",pCutoff=0.05,FCcutoff = 0,col = Volc.cols,cutoffLineType = "blank",legendLabels = Volc.labels))
    
    write.csv(markers,file = paste0(OutPath,'/','Cluster_',unique(extracted.cells.split.cluster[[i]]$seurat_clusters),"_HIVhighpos_vs_HIVneg_Acute_Uninf_dsRed.csv"))
  }
}


#Save of the R workspace#####
print("Saving .RData")
save.image("MAST-DE-FastMNN_11_2023.RData")
print("Done")

#Finish pdf###
dev.off()

##END OF SCRIPT####



