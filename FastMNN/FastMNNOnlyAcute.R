
#Setting the seed to allow reproducibility in terms of UMAPs
set.seed(10403)

#Load necessary packages
library(SeuratWrappers)
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(scCustomize)
library(ggpubr)

#Load the .rds file containing the object.list of seurat objects after the main workflow
#load("/sc/arion/scratch/schmig04/Seurat_BC_test/.RData")

#Acute/Treated/Untreated Merge
#seurat.list.act.tr.untr<-Merge_Seurat_List(seurat.list)

#SCT normalize the merged SeuratObject
#seurat.list.act.tr.untr<-SCTransform(seurat.list.act.tr.untr)

seurat.list.act.tr.untr.MNN<-readRDS("Big_merged_seurat.rds")

#Subset out Uninfected/ART treated

seurat.list.act.MNN <-subset(seurat.list.act.tr.untr.MNN,subset = Condition %in% c("Acute"))

#SCT normalize the merged SeuratObject
seurat.list.act.MNN<-SCTransform(seurat.list.act.MNN)

#Defining variables needed for UMAP
featuresUCOERG<- c("dsRed", "EGFP", "WPRE" )
featuresHIV<- c("gag-pol","pol","pol-vif-vpr-tat-rev","vpu-env","env","mirfp670nano","p2a-cre-ires-nef" )

#Defining variables needed for plotting
HIV.cols <- c('HIV-'='lightgrey','HIV+ low'="#ffd8b1",'HIV+ high'= "purple",'HIV+ very high'="purple")
fluor.cols <- c('Unmarked'='lightgrey','Mixed'='orange','dsRed'='#f5eaea','GFP'='darkgreen')
cond.cols <- c("Acute" = "#cb3e71","ART treated" = "#243684", "Uninfected" ="grey")
mouse.cols <- c("T1" = "#9ebaf5",
                "T2" = "#6d8cd4",
                "T3" = "#1d2e68",
                "A1" = "#ff1919",
                "A2" = "#ff4d4d",
                "A3" = "#ff8080",
                "A0" = "pink",
                "U1" = "#7c7c7c",
                "U2" = "#bcbcbc")

cell.cols<-c("CD4 TCM" = "#F7312E",
             "CD4 Proliferating" = "#0b4aaf",
             "Treg" = "#66A03D",
             "CD4 TEM" = "#CB9801" ,
             "dnT" = "#00A878",
             "CD4 CTL" ="#E06C00",
             "CD4 Naive" = "#88498F",
             "gdT" = "#F5E5FC",
             "CD4 Others" = "#DDE3E3")

cell.order<-c("CD4 Others",
              "gdT",
              "CD4 Naive",
              "CD4 Proliferating",
              "dnT",
              "Treg",
              "CD4 TEM",
              "CD4 CTL",
              "CD4 TCM") #In order of colors


cell.order.rev<-rev(cell.order)

#Perform FastMNN####
#test.features <- VariableFeatures(seurat.list)
#test.features <- setdiff(test.features,featuresHIV)
#mouse.genes<-str_subset(test.features,pattern = "mm10")
#test.features <- setdiff(test.features,mouse.genes)
#test.features <- setdiff(test.features,featuresUCOERG)

#Now we can run the function adn subsequent UMAP

seurat.list.act.MNN <- FindNeighbors(seurat.list.act.MNN,dims = 1:30,reduction = "mnn") 
seurat.list.act.MNN <- FindClusters(seurat.list.act.MNN,resolution = 0.5) 
seurat.list.act.MNN <- RunUMAP(seurat.list.act.MNN,dims = 1:30,reduction = "mnn")


cell.order.abundance<-seurat.list.act.MNN@meta.data%>%
  dplyr::count(predicted.celltype.l2)%>%
  mutate(percent=(n/sum(n))*100)

cell.order.abundance<-arrange(cell.order.abundance,by_group=percent)

cell.order.abundance<-cell.order.abundance$predicted.celltype.l2





pdf("Heatmap2.pdf")

#Seurat Clusters
DimPlot(seurat.list.act.MNN,pt.size = 0.5,reduction = "umap")

DimPlot(seurat.list.act.MNN,pt.size = 0.5)+
  NoLegend()

#Condition
DimPlot(seurat.list.act.MNN,group.by = "Condition",order = c("Acute","ART treated","Uninfected"),pt.size = 0.5)+
  scale_color_manual(values = cond.cols,breaks = c("Acute","ART treated","Uninfected"))

DimPlot(seurat.list.act.MNN,group.by = "Condition",order = c("Acute","ART treated","Uninfected"),pt.size = 0.5)+
  scale_color_manual(values = cond.cols,breaks = c("Acute","ART treated","Uninfected"))+NoLegend()


#MouseID
DimPlot(seurat.list.act.MNN,group.by = "MouseID",cols = mouse.cols,pt.size = 0.5)
DimPlot(seurat.list.act.MNN,group.by = "MouseID",cols = mouse.cols,pt.size = 0.5)+NoLegend()

#HIV expression status
DimPlot(seurat.list.act.MNN,group.by = "status",order = c("HIV+ high","HIV+ low","HIV-"),pt.size = 0.5)+
  scale_color_manual(values = HIV.cols,breaks = c("HIV+ high","HIV+ low","HIV-"))
DimPlot(seurat.list.act.MNN,group.by = "status",order = c("HIV+ high","HIV+ low","HIV-"),pt.size = 0.5)+
  scale_color_manual(values = HIV.cols,breaks = c("HIV+ high","HIV+ low","HIV-"))+NoLegend()

#Predicted celltype
DimPlot(seurat.list.act.MNN,group.by = "predicted.celltype.l2",order = cell.order.abundance,pt.size = 0.5)+
  scale_color_manual(values = cell.cols,breaks = cell.order.rev)

DimPlot(seurat.list.act.MNN,group.by = "predicted.celltype.l2",order = cell.order.abundance,pt.size = 0.5)+
  scale_color_manual(values = cell.cols,breaks = cell.order.rev)

#Fluorescence
DimPlot(seurat.list.act.MNN,group.by = "Fluorescence",pt.size = 0.5,cols = fluor.cols,order = c("GFP","dsRed","Mixed","Unmarked"))

DimPlot(seurat.list.act.MNN,group.by = "Fluorescence",pt.size = 0.5,cols = fluor.cols,order = c("GFP","dsRed","Mixed","Unmarked"))+
  NoLegend()

#Doublet/Singlet
DimPlot(seurat.list.act.MNN,group.by = "Doublet.Singlet",pt.size = 0.5,order = c("Doublet","Singlet"),cols = c("Doublet"="darkblue","Singlet"="lightgrey"))
DimPlot(seurat.list.act.MNN,group.by = "Doublet.Singlet",pt.size = 0.5,order = c("Doublet","Singlet"),cols = c("Doublet"="darkblue","Singlet"="lightgrey"))+
  NoLegend()


#Heatmaps after subsetting
seurat.list.act.MNN<-PrepSCTFindMarkers(seurat.list.act.MNN)

seurat.list.act.MNN.markers <- FindAllMarkers(seurat.list.act.MNN)
mousemarkers<-str_subset(seurat.list.act.MNN.markers$gene,pattern = "mm10")
seurat.list.act.MNN.markers<-seurat.list.act.MNN.markers[!seurat.list.act.MNN.markers$gene %in% mousemarkers,]
seurat.list.act.MNN.markers<-seurat.list.act.MNN.markers[!seurat.list.act.MNN.markers$gene %in% featuresHIV,]
seurat.list.act.MNN.markers<-seurat.list.act.MNN.markers[!seurat.list.act.MNN.markers$gene %in% featuresUCOERG,]

#Markers Top10
seurat.list.act.MNN.markers.top10 <- seurat.list.act.MNN.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

#Show Heatmap
show(DoHeatmap(seurat.list.act.MNN, features = seurat.list.act.MNN.markers.top10$gene,size = 3)+NoLegend()+
       labs(title = "FastMNN-analysis-Heatmap")+theme(text = element_text(size = 5)))


#Various Barplots
h1<-ggplot(seurat.list.act.MNN@meta.data,aes(x=seurat_clusters, fill=predicted.celltype.l2)) + geom_bar(position="fill")+
  scale_y_continuous(labels = scales::percent,expand=c(0,0))+
  theme_bw()+
  labs(fill="Predicted celltype",y="% of cells",x="Cluster",caption="FastMNN - final barplots")+
  theme(axis.text.x = element_text(hjust=1,angle=45),legend.key.size = unit(0.1, 'cm'))+
  scale_fill_manual(values = cell.cols,breaks = cell.order.rev)

show(ggarrange(h1,ncol=1,nrow = 2))


h1<-ggplot(seurat.list.act.MNN@meta.data,aes(x=predicted.celltype.l2, fill=status)) + geom_bar(position="fill")+
  scale_y_continuous(labels = scales::percent,expand=c(0,0))+
  theme_bw()+
  labs(fill="HIV status",y="% of cells",x="predicted.celltypes.l2",caption="FastMNN - final barplots")+
  theme(axis.text.x = element_text(hjust=1,angle=45),legend.key.size = unit(0.1, 'cm'))+
  scale_fill_manual(values = HIV.cols,breaks = c("HIV+ high","HIV+ low","HIV-"))+
  coord_cartesian(ylim = c(0,0.10))

h2<-ggplot(seurat.list.act.MNN@meta.data,aes(x=seurat_clusters, fill=status)) + geom_bar(position="fill")+
  scale_y_continuous(labels = scales::percent,expand=c(0,0))+
  theme_bw()+
  labs(fill="HIV status",y="% of cells",x="predicted.celltypes.l2",caption="FastMNN - final barplots")+
  theme(axis.text.x = element_text(hjust=1,angle=45),legend.key.size = unit(0.1, 'cm'))+
  scale_fill_manual(values = HIV.cols,breaks = c("HIV+ high","HIV+ low","HIV-"))+
  coord_cartesian(ylim = c(0,0.10))
show(ggarrange(h1,h2,ncol=1,nrow = 2))

h1<-ggplot(seurat.list.act.MNN@meta.data,aes(x=predicted.celltype.l2, fill=Condition)) + geom_bar(position="fill")+
  scale_y_continuous(labels = scales::percent,expand=c(0,0))+
  theme_bw()+
  labs(fill="Condition",y="% of cells",x="predicted.celltypes.l2",caption="FastMNN - final barplots")+
  theme(axis.text.x = element_text(hjust=1,angle=45),legend.key.size = unit(0.1, 'cm'))+
  scale_fill_manual(values = cond.cols)

h2<-ggplot(seurat.list.act.MNN@meta.data,aes(x=seurat_clusters, fill=Condition)) + geom_bar(position="fill")+
  scale_y_continuous(labels = scales::percent,expand=c(0,0))+
  theme_bw()+
  labs(fill="Condition",y="% of cells",x="Cluster",caption="FastMNN - final barplots")+
  theme(axis.text.x = element_text(hjust=1,angle=45),legend.key.size = unit(0.1, 'cm'))+
  scale_fill_manual(values = cond.cols)
show(ggarrange(h1,h2,ncol=1,nrow = 2))

o1<-ggplot(seurat.list.act.MNN@meta.data,aes(x=seurat_clusters, fill=Fluorescence))+geom_bar(position="fill")+
  scale_y_continuous(labels = scales::percent,expand=c(0,0))+theme_bw()+
  labs(fill="Fluorescence",y="% of cells",x="Cluster",caption = "FastMNN - final barplots")+
  theme(axis.text.x = element_text(hjust=1,angle=45),legend.key.size = unit(0.1, 'cm'))+
  scale_fill_manual(values = fluor.cols)

o2<-ggplot(seurat.list.act.MNN@meta.data,aes(x=predicted.celltype.l2, fill=Fluorescence)) + geom_bar(position="fill")+
  scale_y_continuous(labels = scales::percent,expand=c(0,0))+theme_bw()+
  labs(fill="Fluorescence",y="% of cells",x="predicted.celltype.l2",caption = "FastMNN - final barplots")+
  theme(axis.text.x = element_text(hjust=1,angle=45),legend.key.size = unit(0.1, 'cm'))+
  scale_fill_manual(values = fluor.cols)

show(ggarrange(o1,o2,ncol=1,nrow = 2))

#
saveRDS(seurat.list.act.MNN,"seurat.list.act.MNN3.rds")

#Subset out Uninfected/ART treated





















#End of script
