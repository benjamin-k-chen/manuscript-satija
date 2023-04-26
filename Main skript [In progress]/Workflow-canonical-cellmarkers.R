#Seurat downstream analysis - WORKFLOW PreSubsetting - 2023-03
#
#
#This is part of the basic workflow for the analysis of the paper Satija et al. (2023)
#
#This workflow was performed on the SeuratObjects pre-subsetting looking at the expression of canonical cellmarker genes
#
#Setting the seed to allow reproducibility in terms of UMAPs
set.seed(10403)

#Load packages we need####
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(DoubletFinder)
library(clustree)
library(scCustomize)

#Load the .rds file
seurat.list <-readRDS("seurat.list.pre.subset.rds")


#Creating the plots
#The list of individual markers can be found in Canonical-cellmakers.R
for (i in 1:length(seurat.list)){
  
  #UMAP DimPlots
  ###
  p1<-DimPlot(seurat.list[[i]])+labs(caption = paste0(dplyr::first(unique(seurat.list[[i]]$Count))," - 1st UMAP"))
  p3<-DimPlot(seurat.list[[i]],group.by = "predicted.celltype.l1",label = TRUE,repel = 10,label.size = 3)+
    labs(caption = paste0(dplyr::first(unique(seurat.list[[i]]$Count))," - 1st UMAP"))+NoLegend()
  p4<-DimPlot(seurat.list[[i]],group.by = "predicted.celltype.l2",label = TRUE,repel = 10,label.size = 3)+
    labs(caption = paste0(dplyr::first(unique(seurat.list[[i]]$Count))," - 1st UMAP"))+NoLegend()
  show(ggarrange(p1,p3,p4,ncol=2,nrow = 2))
  
  
  #VlnPlots
  ##Canonical cell marker features
  #By seurat_clusters
  show(VlnPlot(seurat.list[[i]],features = cellmarker.lymph,group.by = "seurat_clusters",ncol = 2)+
         labs(caption = paste0(dplyr::first(unique(seurat.list[[i]]$Count))," - Lymphocyte markers")))
  
  #By predicted.celltype.l2
  show(VlnPlot(seurat.list[[i]],features = cellmarker.lymph,group.by = "predicted.celltype.l2",ncol = 2)+
         labs(caption = paste0(dplyr::first(unique(seurat.list[[i]]$Count))," - Lymphocyte markers")))
  
  #By seurat_clusters
  show(VlnPlot(seurat.list[[i]],features = cellmarker.immune,group.by = "seurat_clusters",ncol = 2)+
         labs(caption = paste0(dplyr::first(unique(seurat.list[[i]]$Count))," - Immune markers")))
  
  #By predicted.celltype.l2
  show(VlnPlot(seurat.list[[i]],features = cellmarker.immune,group.by = "predicted.celltype.l2",ncol = 2)+
         labs(caption = paste0(dplyr::first(unique(seurat.list[[i]]$Count))," - Immune markers")))
  
  #By seurat_clusters
  show(VlnPlot(seurat.list[[i]],features = cellmarker.cd8,group.by = "seurat_clusters",ncol = 2)+
         labs(caption = paste0(dplyr::first(unique(seurat.list[[i]]$Count))," - CD8 Tcell markers")))
  
  #By predicted.celltype.l2
  show(VlnPlot(seurat.list[[i]],features = cellmarker.cd8,group.by = "predicted.celltype.l2",ncol = 2)+
         labs(caption = paste0(dplyr::first(unique(seurat.list[[i]]$Count))," - CD8 Tcell markers")))
  
  #By seurat_clusters
  show(VlnPlot(seurat.list[[i]],features = cellmarker.bcell,group.by = "seurat_clusters",ncol = 2)+
         labs(caption = paste0(dplyr::first(unique(seurat.list[[i]]$Count))," - Bcell markers")))
  
  #By predicted.celltype.l2
  show(VlnPlot(seurat.list[[i]],features = cellmarker.bcell,group.by = "predicted.celltype.l2",ncol = 2)+
         labs(caption = paste0(dplyr::first(unique(seurat.list[[i]]$Count))," - Bcell markers")))
  
  #By seurat_clusters
  show(VlnPlot(seurat.list[[i]],features = cellmarker.macro,group.by = "seurat_clusters",ncol = 2)+
         labs(caption = paste0(dplyr::first(unique(seurat.list[[i]]$Count))," - Macrophage markers")))
  
  #By predicted.celltype.l2
  show(VlnPlot(seurat.list[[i]],features = cellmarker.macro,group.by = "predicted.celltype.l2",ncol = 2)+
         labs(caption = paste0(dplyr::first(unique(seurat.list[[i]]$Count))," - Macrophage markers")))
  
  #By seurat_clusters
  show(VlnPlot(seurat.list[[i]],features = cellmarker.hspc,group.by = "seurat_clusters",ncol = 2)+
         labs(caption = paste0(dplyr::first(unique(seurat.list[[i]]$Count))," - HSPC markers")))
  
  #By predicted.celltype.l2
  show(VlnPlot(seurat.list[[i]],features = cellmarker.hspc,group.by = "predicted.celltype.l2",ncol = 2)+
         labs(caption = paste0(dplyr::first(unique(seurat.list[[i]]$Count))," - HSPC markers")))
  
  #By seurat_clusters
  show(VlnPlot(seurat.list[[i]],features = cellmarker.dc,group.by = "seurat_clusters",ncol = 2)+
         labs(caption = paste0(dplyr::first(unique(seurat.list[[i]]$Count))," - DC markers")))
  
  #By predicted.celltype.l2
  show(VlnPlot(seurat.list[[i]],features = cellmarker.dc,group.by = "predicted.celltype.l2",ncol = 2)+
         labs(caption = paste0(dplyr::first(unique(seurat.list[[i]]$Count))," - DC markers")))
  
  #By seurat_clusters
  show(VlnPlot(seurat.list[[i]],features = cellmarker.dnT,group.by = "seurat_clusters",ncol = 2)+
         labs(caption = paste0(dplyr::first(unique(seurat.list[[i]]$Count))," - dnT markers")))
  
  #By predicted.celltype.l2
  show(VlnPlot(seurat.list[[i]],features = cellmarker.dnT,group.by = "predicted.celltype.l2",ncol = 2)+
         labs(caption = paste0(dplyr::first(unique(seurat.list[[i]]$Count))," - dnT markers")))
  
  #By seurat_clusters
  show(VlnPlot(seurat.list[[i]],features = cellmarker.gdT,group.by = "seurat_clusters",ncol = 2)+
         labs(caption = paste0(dplyr::first(unique(seurat.list[[i]]$Count))," - gdT markers")))
  
  #By predicted.celltype.l2
  show(VlnPlot(seurat.list[[i]],features = cellmarker.gdT,group.by = "predicted.celltype.l2",ncol = 2)+
         labs(caption = paste0(dplyr::first(unique(seurat.list[[i]]$Count))," - gdT markers")))
  
  #By seurat_clusters
  show(VlnPlot(seurat.list[[i]],features = cellmarker.eryt,group.by = "seurat_clusters",ncol = 2)+
         labs(caption = paste0(dplyr::first(unique(seurat.list[[i]]$Count))," - Erythrocytes markers")))
  
  #By predicted.celltype.l2
  show(VlnPlot(seurat.list[[i]],features = cellmarker.eryt,group.by = "predicted.celltype.l2",ncol = 2)+
         labs(caption = paste0(dplyr::first(unique(seurat.list[[i]]$Count))," - Erythrocytes markers")))
  
  #By seurat_clusters
  show(VlnPlot(seurat.list[[i]],features = cellmarker.platelet,group.by = "seurat_clusters",ncol = 2)+
         labs(caption = paste0(dplyr::first(unique(seurat.list[[i]]$Count))," - Platelet markers")))
  
  #By predicted.celltype.l2
  show(VlnPlot(seurat.list[[i]],features = cellmarker.platelet,group.by = "predicted.celltype.l2",ncol = 2)+
         labs(caption = paste0(dplyr::first(unique(seurat.list[[i]]$Count))," - Platelet markers")))
  
  #By seurat_clusters
  show(VlnPlot(seurat.list[[i]],features = cellmarker.platelet2,group.by = "seurat_clusters",ncol = 2)+
         labs(caption = paste0(dplyr::first(unique(seurat.list[[i]]$Count))," - Platelet markers - 2")))
  
  #By predicted.celltype.l2
  show(VlnPlot(seurat.list[[i]],features = cellmarker.platelet2,group.by = "predicted.celltype.l2",ncol = 2)+
         labs(caption = paste0(dplyr::first(unique(seurat.list[[i]]$Count))," - Platelet markers - 2")))
  
  #By seurat_clusters
  show(VlnPlot(seurat.list[[i]],features = cellmarker.Tfh,group.by = "seurat_clusters",ncol = 2)+
         labs(caption = paste0(dplyr::first(unique(seurat.list[[i]]$Count))," - Tfh markers")))
  
  #By predicted.celltype.l2
  show(VlnPlot(seurat.list[[i]],features = cellmarker.Tfh,group.by = "predicted.celltype.l2",ncol = 2)+
         labs(caption = paste0(dplyr::first(unique(seurat.list[[i]]$Count))," - Tfh markers")))
  
  #By seurat_clusters
  show(VlnPlot(seurat.list[[i]],features = cellmarker.Tfh2,group.by = "seurat_clusters",ncol = 2)+
         labs(caption = paste0(dplyr::first(unique(seurat.list[[i]]$Count))," - Tfh markers - 2")))
  
  #By predicted.celltype.l2
  show(VlnPlot(seurat.list[[i]],features = cellmarker.Tfh2,group.by = "predicted.celltype.l2",ncol = 2)+
         labs(caption = paste0(dplyr::first(unique(seurat.list[[i]]$Count))," - Tfh markers - 2")))
  
  #By seurat_clusters
  show(VlnPlot(seurat.list[[i]],features = cellmarker.TregActive,group.by = "seurat_clusters",ncol = 2)+
         labs(caption = paste0(dplyr::first(unique(seurat.list[[i]]$Count))," - Treg&Activated markers")))
  
  #By predicted.celltype.l2
  show(VlnPlot(seurat.list[[i]],features = cellmarker.TregActive,group.by = "predicted.celltype.l2",ncol = 2)+
         labs(caption = paste0(dplyr::first(unique(seurat.list[[i]]$Count))," - Treg&Activated markers")))
  
  #By seurat_clusters
  show(VlnPlot(seurat.list[[i]],features = cellmarker.Th17cytoly,group.by = "seurat_clusters",ncol = 2)+
         labs(caption = paste0(dplyr::first(unique(seurat.list[[i]]$Count))," - Th17 markers")))
  
  #By predicted.celltype.l2
  show(VlnPlot(seurat.list[[i]],features = cellmarker.Th17cytoly,group.by = "predicted.celltype.l2",ncol = 2)+
         labs(caption = paste0(dplyr::first(unique(seurat.list[[i]]$Count))," - Th17 markers")))
  
  #By seurat_clusters
  show(VlnPlot(seurat.list[[i]],features = cellmarker.EffectorActive,group.by = "seurat_clusters",ncol = 2)+
         labs(caption = paste0(dplyr::first(unique(seurat.list[[i]]$Count))," - Effector activated markers")))
  
  #By predicted.celltype.l2
  show(VlnPlot(seurat.list[[i]],features = cellmarker.EffectorActive,group.by = "predicted.celltype.l2",ncol = 2)+
         labs(caption = paste0(dplyr::first(unique(seurat.list[[i]]$Count))," - Effector activated markers")))
  
  #By seurat_clusters
  show(VlnPlot(seurat.list[[i]],features = cellmarker.Cytotoxic1,group.by = "seurat_clusters",ncol = 2)+
         labs(caption = paste0(dplyr::first(unique(seurat.list[[i]]$Count))," - Cytotoxic markers - 1")))
  
  #By predicted.celltype.l2
  show(VlnPlot(seurat.list[[i]],features = cellmarker.Cytotoxic1,group.by = "predicted.celltype.l2",ncol = 2)+
         labs(caption = paste0(dplyr::first(unique(seurat.list[[i]]$Count))," - Cytotoxic markers -1")))
  
  #By seurat_clusters
  show(VlnPlot(seurat.list[[i]],features = cellmarker.Cytotoxic2,group.by = "seurat_clusters",ncol = 2)+
         labs(caption = paste0(dplyr::first(unique(seurat.list[[i]]$Count))," - Cytotoxic markers - 2")))
  
  #By predicted.celltype.l2
  show(VlnPlot(seurat.list[[i]],features = cellmarker.Cytotoxic2,group.by = "predicted.celltype.l2",ncol = 2)+
         labs(caption = paste0(dplyr::first(unique(seurat.list[[i]]$Count))," - Cytotoxic markers -2")))
  
  #By seurat_clusters
  show(VlnPlot(seurat.list[[i]],features = cellmarker.survival,group.by = "seurat_clusters",ncol = 2)+
         labs(caption = paste0(dplyr::first(unique(seurat.list[[i]]$Count))," - Survival markers - 1")))
  
  #By predicted.celltype.l2
  show(VlnPlot(seurat.list[[i]],features = cellmarker.survival,group.by = "predicted.celltype.l2",ncol = 2)+
         labs(caption = paste0(dplyr::first(unique(seurat.list[[i]]$Count))," - Survival markers -1")))
  
  #By seurat_clusters
  show(VlnPlot(seurat.list[[i]],features = cellmarker.survival2,group.by = "seurat_clusters",ncol = 2)+
         labs(caption = paste0(dplyr::first(unique(seurat.list[[i]]$Count))," - Survival markers - 2")))
  
  #By predicted.celltype.l2
  show(VlnPlot(seurat.list[[i]],features = cellmarker.survival2,group.by = "predicted.celltype.l2",ncol = 2)+
         labs(caption = paste0(dplyr::first(unique(seurat.list[[i]]$Count))," - Survival markers -2")))
  
  #By seurat_clusters
  show(VlnPlot(seurat.list[[i]],features = cellmarker.TCM,group.by = "seurat_clusters",ncol = 2)+
         labs(caption = paste0(dplyr::first(unique(seurat.list[[i]]$Count))," - TCM markers")))
  
  #By predicted.celltype.l2
  show(VlnPlot(seurat.list[[i]],features = cellmarker.TCM,group.by = "predicted.celltype.l2",ncol = 2)+
         labs(caption = paste0(dplyr::first(unique(seurat.list[[i]]$Count))," - TCM markers")))
  
  #By seurat_clusters
  show(VlnPlot(seurat.list[[i]],features = cellmarker.TRM,group.by = "seurat_clusters",ncol = 2)+
         labs(caption = paste0(dplyr::first(unique(seurat.list[[i]]$Count))," - TRM markers")))
  
  #By predicted.celltype.l2
  show(VlnPlot(seurat.list[[i]],features = cellmarker.TRM,group.by = "predicted.celltype.l2",ncol = 2)+
         labs(caption = paste0(dplyr::first(unique(seurat.list[[i]]$Count))," - TRM markers")))
  
  #By seurat_clusters
  show(VlnPlot(seurat.list[[i]],features = cellmarker.activated,group.by = "seurat_clusters",ncol = 2)+
         labs(caption = paste0(dplyr::first(unique(seurat.list[[i]]$Count))," - Activated CD4 markers - 1")))
  
  #By predicted.celltype.l2
  show(VlnPlot(seurat.list[[i]],features = cellmarker.activated,group.by = "predicted.celltype.l2",ncol = 2)+
         labs(caption = paste0(dplyr::first(unique(seurat.list[[i]]$Count))," - Actiavted CD4 markers - 1")))
  
  #By seurat_clusters
  show(VlnPlot(seurat.list[[i]],features = cellmarker.activated2,group.by = "seurat_clusters",ncol = 2)+
         labs(caption = paste0(dplyr::first(unique(seurat.list[[i]]$Count))," - Activated CD4 markers - 2")))
  
  #By predicted.celltype.l2
  show(VlnPlot(seurat.list[[i]],features = cellmarker.activated2,group.by = "predicted.celltype.l2",ncol = 2)+
         labs(caption = paste0(dplyr::first(unique(seurat.list[[i]]$Count))," - Activated CD4 markers - 2")))
  
  #By seurat_clusters
  show(VlnPlot(seurat.list[[i]],features = cellmarker.proliferating,group.by = "seurat_clusters",ncol = 2)+
         labs(caption = paste0(dplyr::first(unique(seurat.list[[i]]$Count))," - Proliferating CD4 markers")))
  
  #By predicted.celltype.l2
  show(VlnPlot(seurat.list[[i]],features = cellmarker.proliferating,group.by = "predicted.celltype.l2",ncol = 2)+
         labs(caption = paste0(dplyr::first(unique(seurat.list[[i]]$Count))," - Proliferating CD4 markers")))
  
  #By seurat_clusters
  show(VlnPlot(seurat.list[[i]],features = cellmarker.Treg,group.by = "seurat_clusters",ncol = 2)+
         labs(caption = paste0(dplyr::first(unique(seurat.list[[i]]$Count))," - Treg markers")))
  
  #By predicted.celltype.l2
  show(VlnPlot(seurat.list[[i]],features = cellmarker.Treg,group.by = "predicted.celltype.l2",ncol = 2)+
         labs(caption = paste0(dplyr::first(unique(seurat.list[[i]]$Count))," - Treg markers")))
  
  #By seurat_clusters
  show(VlnPlot(seurat.list[[i]],features = cellmarker.Th1andTfh,group.by = "seurat_clusters",ncol = 2)+
         labs(caption = paste0(dplyr::first(unique(seurat.list[[i]]$Count))," - Th1 and Tfh markers")))
  
  #By predicted.celltype.l2
  show(VlnPlot(seurat.list[[i]],features = cellmarker.Th1andTfh,group.by = "predicted.celltype.l2",ncol = 2)+
         labs(caption = paste0(dplyr::first(unique(seurat.list[[i]]$Count))," - Th1 anf Tfh markers")))
  
  #By seurat_clusters
  show(VlnPlot(seurat.list[[i]],features = cellmarker.granzyme,group.by = "seurat_clusters",ncol = 2)+
         labs(caption = paste0(dplyr::first(unique(seurat.list[[i]]$Count))," - Cytotoxic T-cell markers")))
  
  #By predicted.celltype.l2
  show(VlnPlot(seurat.list[[i]],features = cellmarker.granzyme,group.by = "predicted.celltype.l2",ncol = 2)+
         labs(caption = paste0(dplyr::first(unique(seurat.list[[i]]$Count))," - Cytotxic T-cell markers")))
  
}


#Based on expression of these canonical cellmarkers as well as the predicted.celltype we get a better understanding of the celltypes in each cluster
#and can consider removal


#END OF SCRIPT#####