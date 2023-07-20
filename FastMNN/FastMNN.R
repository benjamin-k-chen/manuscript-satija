#This is the script used to perform the FastMNN analysis and create the plots

#Setting the seed to allow reproducibility in terms of UMAPs
set.seed(10403)

#Load necessary packages
library(SeuratWrappers)
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(scCustomize)


#Load the .rds file containing the object.list of seurat objects after the main workflow
seurat.list<-readRDS("PATH")

#Acute/Treated/Untreated Merge
seurat.list.act.tr.untr<-Merge_Seurat_List(seurat.list)

#SCT normalize the merged SeuratObject
seurat.list.act.tr.untr<-SCTransform(seurat.list.act.tr.untr)

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
test.features <- VariableFeatures(seurat.list.act.tr.untr)
test.features <- setdiff(test.features,featuresHIV)
mouse.genes<-str_subset(test.features,pattern = "mm10")
test.features <- setdiff(test.features,mouse.genes)
test.features <- setdiff(test.features,featuresUCOERG)

#Create range of resolutions for Clustree
resolution.range <- seq(from = 0, to = 1, by = 0.1)

#Before we can run FastMNN we need to modify the function however, since this version does not
#function with SCT-normalized datasets

##https://github.com/satijalab/seurat/issues/5329 (The error that RunFastMNN will throw)
#this is how they recommend to fix it
trace('RunFastMNN', edit = T)
#replace line 25-28
# objects.sce <- lapply(X = object.list, FUN = function(x,f) { if (DefaultAssay(x) == 'SCT') { x = subset(x = x, features = f) indx <- match(rownames(x@assays$SCT@counts),rownames(x@assays$SCT@scale.data)) x@assays$SCT@scale.data <- x@assays$SCT@scale.data[indx,] }else{ x = subset(x = x, features = f) } return(as.SingleCellExperiment(x)) }, f = features)

#Now we can run the function adn subsequent UMAP

seurat.list.act.tr.untr.MNN <- RunFastMNN(object.list = SplitObject(seurat.list.act.tr.untr,split.by = "Count"), features=test.features)
seurat.list.act.tr.untr.MNN <- RunUMAP(seurat.list.act.tr.untr.MNN,dims = 1:30,reduction = "mnn")
seurat.list.act.tr.untr.MNN <- FindNeighbors(seurat.list.act.tr.untr.MNN,dims = 1:30)
seurat.list.act.tr.untr.MNN <- FindClusters(seurat.list.act.tr.untr.MNN,resolution = 0.5) 

#We get a table with the predicted celltype abundance to use in the order of plotting in the UMAP

cell.order.abundance<-seurat.list.act.tr.untr.MNN@meta.data%>%
   dplyr::count(predicted.celltype.l2)%>%
   mutate(percent=(n/sum(n))*100)

cell.order.abundance<-arrange(cell.order.abundance,by_group=percent)

cell.order.abundance<-cell.order.abundance$predicted.celltype.l2

##Create output pdf and plot the individual plots
pdf("FastMNN.pdf",width = 14,height = 10)

#Each plot is plotted twice, once with legend and once without

#Condition
DimPlot(seurat.list.act.tr.untr.MNN,group.by = "Condition",order = c("Acute","ART treated","Uninfected"),pt.size = 0.5)+
  scale_color_manual(values = cond.cols,breaks = c("Acute","ART treated","Uninfected"))

DimPlot(seurat.list.act.tr.untr.MNN,group.by = "Condition",order = c("Acute","ART treated","Uninfected"),pt.size = 0.5)+
  scale_color_manual(values = cond.cols,breaks = c("Acute","ART treated","Uninfected"))+NoLegend()

#####
#Condition-Version-Acute
cond.cols.acute <- c("Acute" = "#cb3e71","ART treated" ="white", "Uninfected" ="white")

DimPlot(seurat.list.act.tr.untr.MNN,group.by = "Condition",order = c("Acute","ART treated","Uninfected"),pt.size = 0.5)+
  scale_color_manual(values = cond.cols.acute,breaks = c("Acute"))

DimPlot(seurat.list.act.tr.untr.MNN,group.by = "Condition",order = c("Acute","ART treated","Uninfected"),pt.size = 0.5)+
  scale_color_manual(values = cond.cols.acute,breaks = c("Acute"))+NoLegend()

#Condition-Version-ART-treated
cond.cols.treated <- c("Acute" = "white","ART treated" = "#243684", "Uninfected" ="white")

DimPlot(seurat.list.act.tr.untr.MNN,group.by = "Condition",order = c("ART treated","Acute","Uninfected"),pt.size = 0.5)+
  scale_color_manual(values = cond.cols.treated,breaks = c("ART treated"))

DimPlot(seurat.list.act.tr.untr.MNN,group.by = "Condition",order = c("ART treated","Acute","Uninfected"),pt.size = 0.5)+
  scale_color_manual(values = cond.cols.treated,breaks = c("ART treated"))+NoLegend()

#Condition-Version-Uninfected
cond.cols.uninfected <- c("Acute" = "white","ART treated" = "white", "Uninfected" ="grey")

DimPlot(seurat.list.act.tr.untr.MNN,group.by = "Condition",order = c("Uninfected","ART treated","Acute"),pt.size = 0.5)+
  scale_color_manual(values = cond.cols.uninfected,breaks = c("Uninfected"))

DimPlot(seurat.list.act.tr.untr.MNN,group.by = "Condition",order = c("Uninfected","ART treated","Acute"),pt.size = 0.5)+
  scale_color_manual(values = cond.cols.uninfected,breaks = c("Uninfected"))+NoLegend()
########

#MouseID
DimPlot(seurat.list.act.tr.untr.MNN,group.by = "MouseID",cols = mouse.cols,pt.size = 0.5)

DimPlot(seurat.list.act.tr.untr.MNN,group.by = "MouseID",cols = mouse.cols,pt.size = 0.5)+NoLegend()

######
#MouseID-Version-T1
mouse.cols.alt <- c("T1" = "#9ebaf5",
                    "T2" = "white",
                    "T3" = "white",
                    "A1" = "white",
                    "A2" = "white",
                    "A3" = "white",
                    "A0" = "white",
                    "U1" = "white",
                    "U2" = "white")

DimPlot(seurat.list.act.tr.untr.MNN,group.by = "MouseID",order = c("T1"),pt.size = 0.5)+
  scale_color_manual(values = mouse.cols.alt,breaks = c("T1"))

DimPlot(seurat.list.act.tr.untr.MNN,group.by = "MouseID",order = c("T1"),pt.size = 0.5)+
  scale_color_manual(values = mouse.cols.alt,breaks = c("T1"))+NoLegend()

#MouseID-Version-T2
mouse.cols.alt <- c("T1" = "white",
                    "T2" = "#6d8cd4",
                    "T3" = "white",
                    "A1" = "white",
                    "A2" = "white",
                    "A3" = "white",
                    "A0" = "white",
                    "U1" = "white",
                    "U2" = "white")

DimPlot(seurat.list.act.tr.untr.MNN,group.by = "MouseID",order = c("T2"),pt.size = 0.5)+
  scale_color_manual(values = mouse.cols.alt,breaks = c("T2"))

DimPlot(seurat.list.act.tr.untr.MNN,group.by = "MouseID",order = c("T2"),pt.size = 0.5)+
  scale_color_manual(values = mouse.cols.alt,breaks = c("T2"))+NoLegend()

#MouseID-Version-T3
mouse.cols.alt <- c("T1" = "white",
                    "T2" = "white",
                    "T3" = "#1d2e68",
                    "A1" = "white",
                    "A2" = "white",
                    "A3" = "white",
                    "A0" = "white",
                    "U1" = "white",
                    "U2" = "white")

DimPlot(seurat.list.act.tr.untr.MNN,group.by = "MouseID",order = c("T3"),pt.size = 0.5)+
  scale_color_manual(values = mouse.cols.alt,breaks = c("T3"))

DimPlot(seurat.list.act.tr.untr.MNN,group.by = "MouseID",order = c("T3"),pt.size = 0.5)+
  scale_color_manual(values = mouse.cols.alt,breaks = c("T3"))+NoLegend()

#MouseID-Version-A1
mouse.cols.alt <- c("T1" = "white",
                    "T2" = "white",
                    "T3" = "white",
                    "A1" = "#ff1919",
                    "A2" = "white",
                    "A3" = "white",
                    "A0" = "white",
                    "U1" = "white",
                    "U2" = "white")

DimPlot(seurat.list.act.tr.untr.MNN,group.by = "MouseID",order = c("A1"),pt.size = 0.5)+
  scale_color_manual(values = mouse.cols.alt,breaks = c("A1"))

DimPlot(seurat.list.act.tr.untr.MNN,group.by = "MouseID",order = c("A1"),pt.size = 0.5)+
  scale_color_manual(values = mouse.cols.alt,breaks = c("A1"))+NoLegend()

#MouseID-Version-A2
mouse.cols.alt <- c("T1" = "white",
                    "T2" = "white",
                    "T3" = "white",
                    "A1" = "white",
                    "A2" = "#ff4d4d",
                    "A3" = "white",
                    "A0" = "white",
                    "U1" = "white",
                    "U2" = "white")

DimPlot(seurat.list.act.tr.untr.MNN,group.by = "MouseID",order = c("A2"),pt.size = 0.5)+
  scale_color_manual(values = mouse.cols.alt,breaks = c("A2"))

DimPlot(seurat.list.act.tr.untr.MNN,group.by = "MouseID",order = c("A2"),pt.size = 0.5)+
  scale_color_manual(values = mouse.cols.alt,breaks = c("A2"))+NoLegend()

#MouseID-Version-A3
mouse.cols.alt <- c("T1" = "white",
                    "T2" = "white",
                    "T3" = "white",
                    "A1" = "white",
                    "A2" = "white",
                    "A3" = "#ff8080",
                    "A0" = "white",
                    "U1" = "white",
                    "U2" = "white")

DimPlot(seurat.list.act.tr.untr.MNN,group.by = "MouseID",order = c("A3"),pt.size = 0.5)+
  scale_color_manual(values = mouse.cols.alt,breaks = c("A3"))

DimPlot(seurat.list.act.tr.untr.MNN,group.by = "MouseID",order = c("A3"),pt.size = 0.5)+
  scale_color_manual(values = mouse.cols.alt,breaks = c("A3"))+NoLegend()

#MouseID-Version-U1
mouse.cols.alt <- c("T1" = "white",
                    "T2" = "white",
                    "T3" = "white",
                    "A1" = "white",
                    "A2" = "white",
                    "A3" = "white",
                    "A0" = "white",
                    "U1" = "#7c7c7c",
                    "U2" = "white")

DimPlot(seurat.list.act.tr.untr.MNN,group.by = "MouseID",order = c("U1"),pt.size = 0.5)+
  scale_color_manual(values = mouse.cols.alt,breaks = c("U1"))

DimPlot(seurat.list.act.tr.untr.MNN,group.by = "MouseID",order = c("U1"),pt.size = 0.5)+
  scale_color_manual(values = mouse.cols.alt,breaks = c("U1"))+NoLegend()

#MouseID-Version-U1
mouse.cols.alt <- c("T1" = "white",
                    "T2" = "white",
                    "T3" = "white",
                    "A1" = "white",
                    "A2" = "white",
                    "A3" = "white",
                    "A0" = "white",
                    "U1" = "white",
                    "U2" = "#bcbcbc")

DimPlot(seurat.list.act.tr.untr.MNN,group.by = "MouseID",order = c("U2"),pt.size = 0.5)+
  scale_color_manual(values = mouse.cols.alt,breaks = c("U2"))

DimPlot(seurat.list.act.tr.untr.MNN,group.by = "MouseID",order = c("U2"),pt.size = 0.5)+
  scale_color_manual(values = mouse.cols.alt,breaks = c("U2"))+NoLegend()

######


#HIV expression status
DimPlot(seurat.list.act.tr.untr.MNN,group.by = "status",order = c("HIV+ high","HIV+ low","HIV-"),pt.size = 0.5)+
  scale_color_manual(values = HIV.cols,breaks = c("HIV+ high","HIV+ low","HIV-"))

DimPlot(seurat.list.act.tr.untr.MNN,group.by = "status",order = c("HIV+ high","HIV+ low","HIV-"),pt.size = 0.5)+
  scale_color_manual(values = HIV.cols,breaks = c("HIV+ high","HIV+ low","HIV-"))+NoLegend()

#Predicted celltype
DimPlot(seurat.list.act.tr.untr.MNN,group.by = "predicted.celltype.l2",order = cell.order.abundance,pt.size = 0.5)+
  scale_color_manual(values = cell.cols,breaks = cell.order.rev)

DimPlot(seurat.list.act.tr.untr.MNN,group.by = "predicted.celltype.l2",order = cell.order.abundance,pt.size = 0.5)+
  scale_color_manual(values = cell.cols,breaks = cell.order.rev)

#Fluorescence
DimPlot(seurat.list.act.tr.untr.MNN,group.by = "Fluorescence",pt.size = 0.5,cols = fluor.cols,order = c("GFP","dsRed","Mixed","Unmarked"))

DimPlot(seurat.list.act.tr.untr.MNN,group.by = "Fluorescence",pt.size = 0.5,cols = fluor.cols,order = c("GFP","dsRed","Mixed","Unmarked"))+
  NoLegend()

#Doublet/Singlet
DimPlot(seurat.list.act.tr.untr.MNN,group.by = "Doublet.Singlet",pt.size = 0.5,order = c("Doublet","Singlet"),cols = c("Doublet"="darkblue","Singlet"="lightgrey"))
DimPlot(seurat.list.act.tr.untr.MNN,group.by = "Doublet.Singlet",pt.size = 0.5,order = c("Doublet","Singlet"),cols = c("Doublet"="darkblue","Singlet"="lightgrey"))+
  NoLegend()



##Combined UMI plot

seurat.list.act.tr.untr.MNN.AandT<-subset(seurat.list.act.tr.untr.MNN,subset = MouseID %in% c("A0","A1","A2","A3","T1","T2","T3")) #Only using the treated and acuted mice

gene_HIV2 <- ifelse(str_detect(rownames(seurat.list.act.tr.untr.MNN.AandT@assays$RNA),"gag-pol|^pol|pol-vif-vpr-tat-rev|vpu-env|env|p2a-cre-ires-nef"),"HIV+", "HIV-")
hivpos_inds2 <- gene_HIV2 == "HIV+"
hivneg_inds2 <- gene_HIV2 == "HIV-"
cell_hivstat2 <- tibble(n_hivpos_umi = Matrix::colSums(seurat.list.act.tr.untr.MNN.AandT@assays$RNA[hivpos_inds2,]),
                        n_hivneg_umi = Matrix::colSums(seurat.list.act.tr.untr.MNN.AandT@assays$RNA[hivneg_inds2,]),
                        tot_umi = Matrix::colSums(seurat.list.act.tr.untr.MNN.AandT@assays$RNA),
                        Status = seurat.list.act.tr.untr.MNN.AandT$status,
                        Count = seurat.list.act.tr.untr.MNN.AandT$Count,
                        Condition = seurat.list.act.tr.untr.MNN.AandT$Condition,
                        MouseID = seurat.list.act.tr.untr.MNN.AandT$MouseID,
                        Mean_hivpos_umi = base::mean(n_hivpos_umi))
cell_hivstat2<-cell_hivstat2%>%
  arrange(Status)
ggplot(cell_hivstat2, aes(n_hivpos_umi, n_hivneg_umi, color = Status,shape=MouseID)) +
  geom_point(size=3)+scale_x_log10()+scale_y_log10()+labs(title = "HIV+ vs. HIV- UMIs")+
  scale_color_manual(values = HIV.cols,breaks = c("HIV+ high","HIV+ low","HIV-"))+
  theme_bw(base_size = 20)

#Finishing pdf
dev.off()

##Barplots with HIV+ cell throughout seurat clusters
pdf("Barplots-FastMNN.pdf")
h1<-ggplot(seurat.list.act.tr.untr.MNN@meta.data,aes(x=predicted.celltype.l2, fill=status)) + geom_bar(position="fill")+
  scale_y_continuous(labels = scales::percent,expand=c(0,0))+
  theme_bw()+
  labs(fill="HIV status",y="% of cells",x="predicted.celltypes.l2",caption="FastMNN - final barplots")+
  theme(axis.text.x = element_text(hjust=1,angle=45),legend.key.size = unit(0.1, 'cm'))+
  scale_fill_manual(values = HIV.cols,breaks = c("HIV+ high","HIV+ low","HIV-"))+
  coord_cartesian(ylim = c(0,0.05))

h2<-ggplot(seurat.list.act.tr.untr.MNN@meta.data,aes(x=seurat_clusters, fill=status)) + geom_bar(position="fill")+
  scale_y_continuous(labels = scales::percent,expand=c(0,0))+
  theme_bw()+
  labs(fill="HIV status",y="% of cells",x="predicted.celltypes.l2",caption="FastMNN - final barplots")+
  theme(axis.text.x = element_text(hjust=1,angle=45),legend.key.size = unit(0.1, 'cm'))+
  scale_fill_manual(values = HIV.cols,breaks = c("HIV+ high","HIV+ low","HIV-"))+
  coord_cartesian(ylim = c(0,0.10))
show(ggarrange(h1,h2,ncol=1,nrow = 2))

h1<-ggplot(seurat.list.act.tr.untr.MNN@meta.data,aes(x=predicted.celltype.l2, fill=Condition)) + geom_bar(position="fill")+
  scale_y_continuous(labels = scales::percent,expand=c(0,0))+
  theme_bw()+
  labs(fill="Condition",y="% of cells",x="predicted.celltypes.l2",caption="FastMNN - final barplots")+
  theme(axis.text.x = element_text(hjust=1,angle=45),legend.key.size = unit(0.1, 'cm'))+
  scale_fill_manual(values = cond.cols)

h2<-ggplot(seurat.list.act.tr.untr.MNN@meta.data,aes(x=seurat_clusters, fill=Condition)) + geom_bar(position="fill")+
  scale_y_continuous(labels = scales::percent,expand=c(0,0))+
  theme_bw()+
  labs(fill="Condition",y="% of cells",x="Cluster",caption="FastMNN - final barplots")+
  theme(axis.text.x = element_text(hjust=1,angle=45),legend.key.size = unit(0.1, 'cm'))+
  scale_fill_manual(values = cond.cols)
show(ggarrange(h1,h2,ncol=1,nrow = 2))

o1<-ggplot(seurat.list.act.tr.untr.MNN@meta.data,aes(x=seurat_clusters, fill=Fluorescence))+geom_bar(position="fill")+
  scale_y_continuous(labels = scales::percent,expand=c(0,0))+theme_bw()+
  labs(fill="Fluorescence",y="% of cells",x="Cluster",caption = "FastMNN - final barplots")+
  theme(axis.text.x = element_text(hjust=1,angle=45),legend.key.size = unit(0.1, 'cm'))+
  scale_fill_manual(values = fluor.cols)

o2<-ggplot(seurat.list.act.tr.untr.MNN@meta.data,aes(x=predicted.celltype.l2, fill=Fluorescence)) + geom_bar(position="fill")+
  scale_y_continuous(labels = scales::percent,expand=c(0,0))+theme_bw()+
  labs(fill="Fluorescence",y="% of cells",x="predicted.celltype.l2",caption = "FastMNN - final barplots")+
  theme(axis.text.x = element_text(hjust=1,angle=45),legend.key.size = unit(0.1, 'cm'))+
  scale_fill_manual(values = fluor.cols)

show(ggarrange(o1,o2,ncol=1,nrow = 2))
dev.off()


#END OF SCRIPT####
