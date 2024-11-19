
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

#Acute/Treated/Untreated Merge
seurat.list.act.tr.untr<-Merge_Seurat_List(seurat.list)

#SCT normalize the merged SeuratObject
seurat.list.act.tr.untr<-SCTransform(seurat.list.act.tr.untr)

seurat.list.act.tr.untr.MNN<-readRDS("Big_merged_seurat.rds")

#Subset out Uninfected/ART treated
seurat.list.act.MNN <-subset(seurat.list.act.tr.untr.MNN,subset = Condition %in% c("Acute"))

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

pdf("Acute.pdf")

#Seurat Clusters
DimPlot(seurat.list.act.MNN,pt.size = 0.5,order = rev(cell.orders.clusters))
DimPlot(seurat.list.act.MNN,pt.size = 0.5,order = rev(cell.orders.clusters))+NoLegend()

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
saveRDS(seurat.list.act.MNN,"seurat.list.act.MNN.rds")


##############

#Let's extract the numbers of cells HIV-/HIV+ using cluster designations

total_vector<-seurat.list.act.MNN@meta.data%>%
  group_by(seurat_clusters)%>%
  dplyr::count(status=="HIV-")

HIVpos.table<-total_vector[total_vector$`status == "HIV-"`=="FALSE",]
HIVneg.table<-total_vector[total_vector$`status == "HIV-"`=="TRUE",]

total_vector1<-seurat.list.act.MNN@meta.data%>%
  group_by(seurat_clusters)%>%
  dplyr::count(seurat_clusters)

TableA <- data.frame(
  HIVpos=c(total_vector1$n-HIVneg.table$n),#Sucesses
  HIVneg=HIVneg.table$n#Failures
)

TableA$RowSum<-rowSums(TableA)

#Unfortunately I'm using a version of the seurat object which doesn't have the final cluster desigantions so I'll add them manually

TableA$Clusters<-0:13


TableA<-`rownames<-`(TableA,c("Th1/Th17 Memory",
                              "Proliferating CD4",
                              "Tfh Memory",
                              "Tregs",
                              "IL10RA+ Tregs",
                              "Th1 Proliferating CD4",
                              "Tfh/Th22",
                              "MT Cells",
                              "Naïve CD4-1",
                              "Naïve CD4-2",
                              "Cytolytic Memory",
                              "Cytolytic Proliferating CD4",
                              "Th2",
                              "Th22/Th9 Memory"
))

TableA$CellType<-rownames(TableA)

TableA_Total<-data_frame(HIVpos=sum(TableA$HIVpos),HIVneg=sum(TableA$HIVneg))

###HIV Barplot
cell.orders<-c("Proliferating CD4",
               "Th1 Proliferating CD4",
               "Cytolytic Proliferating CD4",
               "Naïve CD4-1",
               "Naïve CD4-2",
               "Cytolytic Memory",
               "Tfh Memory",
               "Th1/Th17 Memory",
               "Th22/Th9 Memory",
               "Tfh/Th22",
               "Th2",
               "Tregs",
               "IL10RA+ Tregs",
               "MT Cells")

cell.orders.clusters<-c("1","5","11","8","9","10","2","0","13","6","12","3","4","7")



ggplot(seurat.list.act.MNN@meta.data,aes(x=seurat_clusters, fill=status))+
  geom_bar(position="fill")+
  scale_y_continuous(labels = scales::percent,expand=c(0,0))+
  scale_x_discrete(limits=cell.orders.clusters,labels=cell.orders)+
  theme_bw()+
  labs(fill="HIV status",y="% of cells",x="")+
  theme(axis.text.x = element_text(hjust=1,angle=45),legend.key.size = unit(0.1, 'cm'))+
  coord_cartesian(ylim = c(0,0.1))+
  scale_fill_manual(values = HIV.cols,breaks = c("HIV+ high","HIV+ low","HIV-"))+
  geom_hline(yintercept = TableA_Total$HIVpos/(TableA_Total$HIVneg+TableA_Total$HIVpos), color="black",linetype=2)


dev.off()

##########Treated

seurat.list.tr.MNN <-subset(`seurat.list.act.mtr.untr.MNN.clusters copy`,subset = MouseID %in% c("T1","T2"))

pdf("Treated.pdf")

#Seurat Clusters
DimPlot(seurat.list.tr.MNN,pt.size = 0.5,reduction = "umap",order = rev(cell.orders.clusters))

DimPlot(seurat.list.tr.MNN,pt.size = 0.5,order = rev(cell.orders.clusters))+
  NoLegend()
#No Th2 in Treated

#Condition
DimPlot(seurat.list.tr.MNN,group.by = "Condition",order = c("Acute","ART treated","Uninfected"),pt.size = 0.5)+
  scale_color_manual(values = cond.cols,breaks = c("Acute","ART treated","Uninfected"))

DimPlot(seurat.list.tr.MNN,group.by = "Condition",order = c("Acute","ART treated","Uninfected"),pt.size = 0.5)+
  scale_color_manual(values = cond.cols,breaks = c("Acute","ART treated","Uninfected"))+NoLegend()

#MouseID
DimPlot(seurat.list.tr.MNN,group.by = "MouseID",cols = mouse.cols,pt.size = 0.5)
DimPlot(seurat.list.tr.MNN,group.by = "MouseID",cols = mouse.cols,pt.size = 0.5)+NoLegend()

#HIV expression status
DimPlot(seurat.list.tr.MNN,group.by = "status",order = c("HIV+ high","HIV+ low","HIV-"),pt.size = 0.5)+
  scale_color_manual(values = HIV.cols,breaks = c("HIV+ high","HIV+ low","HIV-"))
DimPlot(seurat.list.tr.MNN,group.by = "status",order = c("HIV+ high","HIV+ low","HIV-"),pt.size = 0.5)+
  scale_color_manual(values = HIV.cols,breaks = c("HIV+ high","HIV+ low","HIV-"))+NoLegend()

#Predicted celltype
DimPlot(seurat.list.tr.MNN,group.by = "predicted.celltype.l2",order = cell.order.abundance,pt.size = 0.5)+
  scale_color_manual(values = cell.cols,breaks = cell.order.rev)

DimPlot(seurat.list.tr.MNN,group.by = "predicted.celltype.l2",order = cell.order.abundance,pt.size = 0.5)+
  scale_color_manual(values = cell.cols,breaks = cell.order.rev)

#Fluorescence
DimPlot(seurat.list.tr.MNN,group.by = "Fluorescence",pt.size = 0.5,cols = fluor.cols,order = c("GFP","dsRed","Mixed","Unmarked"))

DimPlot(seurat.list.tr.MNN,group.by = "Fluorescence",pt.size = 0.5,cols = fluor.cols,order = c("GFP","dsRed","Mixed","Unmarked"))+
  NoLegend()

#Various Barplots
total_vector<-seurat.list.tr.MNN@meta.data%>%
  group_by(seurat_clusters)%>%
  dplyr::count(Fluorescence=="dsRed")

GFP.table<-total_vector[total_vector$`Fluorescence == "dsRed"`=="FALSE",]
dsRed.table<-total_vector[total_vector$`Fluorescence == "dsRed"`=="TRUE",]

total_vector1<-seurat.list.tr.MNN@meta.data%>%
  group_by(seurat_clusters)%>%
  dplyr::count(seurat_clusters)

TableE <- data.frame(
  GFP=c(total_vector1$n-dsRed.table$n),#Sucesses
  dsRed=dsRed.table$n#Failures
)

#Unfortunately I'm using a version of the seurat object which doesn't have the final cluster desigantions so I'll add them manually
TableE<-`rownames<-`(TableE,c("Th1/Th17 Memory",
                              "Proliferating CD4", 
                              "Tfh Memory",
                                                          "Tregs",
                                                          "IL10RA+ Tregs",
                                                          "Th1 Proliferating CD4",
                                                          "Tfh/Th22",
                                                          "MT Cells",
                                                          "Naïve CD4-1",
                                                          "Naïve CD4-2",
                                                          "Cytolytic Memory",
                                                          "Cytolytic Proliferating CD4",
                                                          "Th22/Th9 Memory"
))


TableE_Total<-data_frame(GFP=sum(TableE$GFP),dsRed=sum(TableE$dsRed))

ggplot(seurat.list.tr.MNN@meta.data,aes(x=seurat_clusters, fill=Fluorescence))+geom_bar(position="fill")+
  scale_y_continuous(labels = scales::percent,expand=c(0,0))+theme_bw()+
  scale_x_discrete(limits=cell.orders.clusters,labels=cell.orders)+
  coord_cartesian(ylim = c(0,0.05))+
  labs(fill="Fluorescence",y="% of cells",x="Cluster",caption = "FastMNN - final barplots")+
  theme(axis.text.x = element_text(hjust=1,angle=45),legend.key.size = unit(0.1, 'cm'))+
  scale_fill_manual(values = fluor.cols)+
  geom_hline(yintercept = TableE_Total$GFP/(TableE_Total$dsRed+TableE_Total$GFP), color="black",linetype=2)

TableE_fisher<-fisher_test(TableE,simulate.p.value = TRUE)

b1<-ggplot(TableA, aes(fill=rownames(TableA), y=rowSums(TableA), x="")) + 
  geom_bar(position="fill", stat="identity")+
  coord_polar("y", start=0)+
  scale_y_continuous(labels = scales::percent,expand=c(0,0))+
  labs(fill="",y="% of cells",x="Acute")+
  theme_void()

b2<-ggplot(TableE, aes(fill=rownames(TableE), y=rowSums(TableE), x="")) + 
  geom_bar(position="fill", stat="identity")+
  coord_polar("y", start=0)+
  scale_y_continuous(labels = scales::percent,expand=c(0,0))+
  labs(fill="",y="% of cells",x="ART treated")+
  theme_void()

ggarrange(b1,b2)

pdf("BarplotPopulations.pdf")
b3<-ggplot(TableA, aes(fill=CellType,y=RowSum, x="")) + 
  geom_bar(position = "fill",stat = "identity")+
  scale_fill_discrete(labels=cell.orders,breaks=cell.orders,labels=CellType)+
  scale_y_continuous(labels = scales::percent,expand=c(0,0))+
  labs(fill="",y="% of cells",x="Acute")+
  theme_bw()

b4<-ggplot(TableE2, aes(fill=rownames(TableE2), y=rowSums(TableE2), x="")) + 
  geom_bar(position="fill", stat="identity")+
  scale_y_continuous(labels = scales::percent,expand=c(0,0))+
  labs(fill="",y="% of cells",x="ART treated")+
  theme_bw()

ggarrange(b3,b4,ncol=2,nrow = 1)

dev.off()

#End of script