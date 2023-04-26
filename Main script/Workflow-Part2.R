#Seurat downstream analysis - BASIC WORKFLOW - Part 2 - Subsetting and final plots
#
#
#This is the basic workflow for the analysis of the paper Satija et al. (2023)
#
#
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

#Loading the dataset####
#Load the .rds file
seurat.list <-readRDS("seurat.list.pre.subset.rds")


#Subset the datasets - based on plots from Workflow-canonical cellmarkers
`%notin%` <- Negate(`%in%`)

#Subsetting
for (i in 1:length(seurat.list)){
  print("Subsetting datasets")
  if  (dplyr::first(seurat.list[[i]]$Count) %in% c("27M1")){seurat.list[[i]]<- subset(seurat.list[[i]], subset = seurat_clusters %in% c(2,3,4))}
  
  if  (dplyr::first(seurat.list[[i]]$Count) %in% c("297_count")){seurat.list[[i]]<- subset(seurat.list[[i]], subset = predicted.celltype.l2 %notin% c("CD8 TCM","CD8 TEM","B intermediate"))}
  
  if  (dplyr::first(seurat.list[[i]]$Count) %in% c("368_count")){seurat.list[[i]]<- subset(seurat.list[[i]], subset = seurat_clusters %notin% c(10)
                                                                                                                                                   &predicted.celltype.l2 %notin% c("CD16 Mono"))}
  
  if  (dplyr::first(seurat.list[[i]]$Count) %in% c("609_count")){seurat.list[[i]]<- subset(seurat.list[[i]], subset=seurat_clusters %notin% c(6,9)
                                                                                                                                                   &predicted.celltype.l2 %notin% c("HSPC"))}
  
  if  (dplyr::first(seurat.list[[i]]$Count) %in% c("611_count")){seurat.list[[i]]<- subset(seurat.list[[i]], subset=seurat_clusters %notin% c(7,6)
                                                                                                                                                   &predicted.celltype.l2 %notin% c("CD14 Mono","CD16 Mono"))}
  
  if  (dplyr::first(seurat.list[[i]]$Count) %in% c("618_count")){seurat.list[[i]]<- subset(seurat.list[[i]], subset=seurat_clusters %notin% c(8))}
  
  if  (dplyr::first(seurat.list[[i]]$Count) %in% c("620_count")){seurat.list[[i]]<- subset(seurat.list[[i]], subset=seurat_clusters %notin% c(4,9))}
  
}

#Normalization
seurat.list<-lapply(X=seurat.list,FUN = SCTransform,vars.to.regress = "mitoHuCH38Percent")


#We confirmed several misnomers in the datasets and renamed them as 'CD4 Others')
#Renaming predicted.celltypes.l2
print("Consolidating predicted.celltype.l2")

for(i in 1:length(seurat.list)) {
  seurat.list[[i]]$predicted.celltype.l2old<-seurat.list[[i]]$predicted.celltype.l2
  seurat.list[[i]]$predicted.celltype.l2 <- replace(seurat.list[[i]]$predicted.celltype.l2,c("CD8 Naive|CD8 Proliferating|CD8 TEM||CD8 TCM|CD14 Mono|CD16 Mono|cDC1|cDC2|pDC|B intermediate|B naive|MAIT|HSPC|Platelet|B naive|B intermediate|NK|Plasmablast|Eryth|NK Proliferating|ILC"),"CD4 Others")
  seurat.list[[i]]$predicted.celltype.l2 <- replace(seurat.list[[i]]$predicted.celltype.l2,"CD4 Others Proliferating^","CD4 Others")
}


#Final UMAPs
datasetsUMAPtbl<-data.frame(datasetsUMAP2)
datasetsUMAPtbl$Count<-rownames(datasetsUMAPtbl)

datasetsClusterRestbl<-data.frame(datasetsClusterRes2)
datasetsClusterRestbl$Count<-rownames(datasetsClusterRestbl)

datasetsUMAPandClusterRestbl<-inner_join(datasetsUMAPtbl,datasetsClusterRestbl,by="Count")



for (i in 1:length(seurat.list)){
    
    tt2sub<-seurat.list[[i]]@meta.data%>%
      dplyr::count(predicted.celltype.l2)%>%
      mutate(percent=(n/sum(n))*100)
    
    tt2sub2<-arrange(tt2sub,by_group=percent)
    
    
    seurat.info <- tibble(Count = seurat.list[[i]]$Count)
    seurat.info <- inner_join(seurat.info,datasetsUMAPandClusterRestbl,by="Count")
    
    test.features <- VariableFeatures(seurat.list[[i]])
    test.features <- setdiff(test.features,featuresHIV)
    mouse.genes<-str_subset(test.features,pattern = "mm10")
    test.features <- setdiff(test.features,mouse.genes)
    test.features <- setdiff(test.features,featuresUCOERG)
    seurat.list[[i]] <- RunPCA(seurat.list[[i]],features=test.features)
    seurat.list[[i]] <- RunUMAP(seurat.list[[i]],dims = 1:unique(seurat.info$datasetsUMAP2))
    seurat.list[[i]] <- RunTSNE(seurat.list[[i]],dims = 1:unique(seurat.info$datasetsUMAP2))
    seurat.list[[i]] <- FindNeighbors(seurat.list[[i]],dims = 1:unique(seurat.info$datasetsUMAP2))
    seurat.list[[i]] <- FindClusters(seurat.list[[i]],resolution = unique(seurat.info$datasetsClusterRes2)) 
    
    p1<-DimPlot(seurat.list[[i]])+NoLegend()+labs(title = paste0(dplyr::first(seurat.list[[i]]$Count),"-UMAP - post-subsetting"),subtitle = paste0("ndims=1:",dplyr::first(seurat.info$datasetsUMAP2)," res=",dplyr::first(seurat.info$datasetsClusterRes2)))
    p2<-DimPlot(seurat.list[[i]],group.by = "Count")+NoLegend()+labs(title = paste0(dplyr::first(seurat.list[[i]]$Count),"-UMAP - post-subsetting"),subtitle = paste0("ndims=1:",dplyr::first(seurat.info$datasetsUMAP2)," res=",dplyr::first(seurat.info$datasetsClusterRes2)))
    p3<-DimPlot(seurat.list[[i]],group.by = "predicted.celltype.l1")+NoLegend()+labs(title = paste0(dplyr::first(seurat.list[[i]]$Count),"-UMAP - post-subsetting"),subtitle = paste0("ndims=1:",dplyr::first(seurat.info$datasetsUMAP2)," res=",dplyr::first(seurat.info$datasetsClusterRes2)))
    p4<-DimPlot(seurat.list[[i]],group.by = "predicted.celltype.l2",order=tt2sub2$predicted.celltype.l2)+NoLegend()+labs(title = paste0(dplyr::first(seurat.list[[i]]$Count),"-UMAP - post-subsetting"),subtitle = paste0("ndims=1:",dplyr::first(seurat.info$datasetsUMAP2)," res=",dplyr::first(seurat.info$datasetsClusterRes2)))+
      scale_color_manual(values = cell.cols,breaks = cell.order.rev)
    show(ggarrange(p1,p2,p3,p4,ncol=2,nrow = 2))
    
    p1<-DimPlot(seurat.list[[i]])+labs(title = paste0(dplyr::first(seurat.list[[i]]$Count),"-UMAP - post-subsetting"),subtitle = paste0("ndims=1:",dplyr::first(seurat.info$datasetsUMAP2)," res=",dplyr::first(seurat.info$datasetsClusterRes2)))
    p2<-DimPlot(seurat.list[[i]],group.by = "Count")+labs(title = paste0(dplyr::first(seurat.list[[i]]$Count),"-UMAP - post-subsetting"),subtitle = paste0("ndims=1:",dplyr::first(seurat.info$datasetsUMAP2)," res=",dplyr::first(seurat.info$datasetsClusterRes2)))
    p3<-DimPlot(seurat.list[[i]],group.by = "predicted.celltype.l1")+labs(title = paste0(dplyr::first(seurat.list[[i]]$Count),"-UMAP - post-subsetting"),subtitle = paste0("ndims=1:",dplyr::first(seurat.info$datasetsUMAP2)," res=",dplyr::first(seurat.info$datasetsClusterRes2)))
    p4<-DimPlot(seurat.list[[i]],group.by = "predicted.celltype.l2",order=tt2sub2$predicted.celltype.l2)+labs(title = paste0(dplyr::first(seurat.list[[i]]$Count),"-UMAP - post-subsetting"),subtitle = paste0("ndims=1:",dplyr::first(seurat.info$datasetsUMAP2)," res=",dplyr::first(seurat.info$datasetsClusterRes2)))+
      scale_color_manual(values = cell.cols,breaks = cell.order.rev)
    show(ggarrange(p1,p2,p3,p4,ncol=2,nrow = 2))
    
    p5<-DimPlot(seurat.list[[i]],group.by = "Condition",cols = cond.cols)+NoLegend()+labs(title = paste0(dplyr::first(seurat.list[[i]]$Count),"-UMAP - post-subsetting"),subtitle = paste0("ndims=1:",dplyr::first(seurat.info$datasetsUMAP2)," res=",dplyr::first(seurat.info$datasetsClusterRes2)))
    p6<-DimPlot(seurat.list[[i]],group.by = "status",order=c("HIV+ high","HIV+ low","HIV-"),cols = HIV.cols)+NoLegend()+labs(title = paste0(dplyr::first(seurat.list[[i]]$Count),"-UMAP - post-subsetting"),subtitle = paste0("ndims=1:",dplyr::first(seurat.info$datasetsUMAP2)," res=",dplyr::first(seurat.info$datasetsClusterRes2)))
    p7<-DimPlot(seurat.list[[i]],group.by = "Fluorescence",order=c("GFP","dsRed","Mixed","Unmarked"),cols = fluor.cols)+NoLegend()+labs(title = paste0(dplyr::first(seurat.list[[i]]$Count),"-UMAP - post-subsetting"),subtitle = paste0("ndims=1:",dplyr::first(seurat.info$datasetsUMAP2)," res=",dplyr::first(seurat.info$datasetsClusterRes2)))
    p8<-DimPlot(seurat.list[[i]],group.by = "LentiRG.expression",order=c("LentiRG","No epxression"),cols = dsRed.cols)+NoLegend()+labs(title = paste0(dplyr::first(seurat.list[[i]]$Count),"-UMAP - post-subsetting"),subtitle = paste0("ndims=1:",dplyr::first(seurat.info$datasetsUMAP2)," res=",dplyr::first(seurat.info$datasetsClusterRes2)))
    show(ggarrange(p6,p7,p8,ncol=2,nrow = 2))
    
    p5<-DimPlot(seurat.list[[i]],group.by = "Condition",cols = cond.cols)+labs(title = paste0(dplyr::first(seurat.list[[i]]$Count),"-UMAP - post-subsetting"),subtitle = paste0("ndims=1:",dplyr::first(seurat.info$datasetsUMAP2)," res=",dplyr::first(seurat.info$datasetsClusterRes2)))
    p6<-DimPlot(seurat.list[[i]],group.by = "status",order=c("HIV+ high","HIV+ low","HIV-"),cols = HIV.cols)+labs(title = paste0(dplyr::first(seurat.list[[i]]$Count),"-UMAP - post-subsetting"),subtitle = paste0("ndims=1:",dplyr::first(seurat.info$datasetsUMAP2)," res=",dplyr::first(seurat.info$datasetsClusterRes2)))
    p7<-DimPlot(seurat.list[[i]],group.by = "Fluorescence",order=c("GFP","dsRed","Mixed","Unmarked"),cols = fluor.cols)+labs(title = paste0(dplyr::first(seurat.list[[i]]$Count),"-UMAP - post-subsetting"),subtitle = paste0("ndims=1:",dplyr::first(seurat.info$datasetsUMAP2)," res=",dplyr::first(seurat.info$datasetsClusterRes2)))
    p8<-DimPlot(seurat.list[[i]],group.by = "LentiRG.expression",order=c("LentiRG","No epxression"),cols = dsRed.cols)+labs(title = paste0(dplyr::first(seurat.list[[i]]$Count),"-UMAP - post-subsetting"),subtitle = paste0("ndims=1:",dplyr::first(seurat.info$datasetsUMAP2)," res=",dplyr::first(seurat.info$datasetsClusterRes2)))
    show(ggarrange(p6,p7,p8,ncol=2,nrow = 2))
    
    
    p9<-DimPlot(seurat.list[[i]],group.by = "MouseID",cols = mouse.cols)+NoLegend()+labs(title = paste0(dplyr::first(seurat.list[[i]]$Count),"-UMAP - post-subsetting"),subtitle = paste0("ndims=1:",dplyr::first(seurat.info$datasetsUMAP2)," res=",dplyr::first(seurat.info$datasetsClusterRes2)))
    p10<-DimPlot(seurat.list[[i]],group.by = "predicted.celltype.l2old")+NoLegend()+labs(title = paste0(dplyr::first(seurat.list[[i]]$Count),"-UMAP - post-subsetting"),subtitle = paste0("ndims=1:",dplyr::first(seurat.info$datasetsUMAP2)," res=",dplyr::first(seurat.info$datasetsClusterRes2)))
    show(ggarrange(p9,p10,ncol=2,nrow = 2))
    
    p9<-DimPlot(seurat.list[[i]],group.by = "MouseID",cols = mouse.cols)+labs(title = paste0(dplyr::first(seurat.list[[i]]$Count),"-UMAP - post-subsetting"),subtitle = paste0("ndims=1:",dplyr::first(seurat.info$datasetsUMAP2)," res=",dplyr::first(seurat.info$datasetsClusterRes2)))
    p10<-DimPlot(seurat.list[[i]],group.by = "predicted.celltype.l2old")+labs(title = paste0(dplyr::first(seurat.list[[i]]$Count),"-UMAP - post-subsetting"),subtitle = paste0("ndims=1:",dplyr::first(seurat.info$datasetsUMAP2)," res=",dplyr::first(seurat.info$datasetsClusterRes2)))
    show(ggarrange(p9,p10,ncol=2,nrow = 2))
    
    #Heatmaps after subsetting
    seurat.list.markers <- FindAllMarkers(seurat.list[[i]])
    mousemarkers<-str_subset(seurat.list.markers$gene,pattern = "mm10")
    seurat.list.markers<-seurat.list.markers[!seurat.list.markers$gene %in% mousemarkers,]
    seurat.list.markers<-seurat.list.markers[!seurat.list.markers$gene %in% featuresHIV,]
    seurat.list.markers<-seurat.list.markers[!seurat.list.markers$gene %in% featuresUCOERG,]
    
    #Save Markers full list and Top10
    assign(paste0(dplyr::first(seurat.list[[i]]$Count),"-markers-unsupervised"),seurat.list.markers)
    seurat.list.markers.top10 <- seurat.list.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
    assign(paste0(dplyr::first(seurat.list[[i]]$Count),"-markers-unsupervised-topten"),seurat.list.markers.top10)
    
    #Show Heatmap
    show(DoHeatmap(seurat.list[[i]], features = seurat.list.markers.top10$gene,size = 3)+NoLegend()+
           labs(title = paste0(dplyr::first(seurat.list[[i]]$Count)),
                subtitle = "After cluster subsetting")+theme(text = element_text(size = 5)))
    
}


#Save the workspace and the .rds file

save.image("Seurat_final.RData")

saveRDS(seurat.list,"seurat_list_final.rds")



#END OF SCRIPT#



