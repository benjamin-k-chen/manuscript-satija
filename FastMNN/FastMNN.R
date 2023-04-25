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

test.features <- VariableFeatures(seurat)
test.features <- setdiff(test.features,featuresHIV)
mouse.genes<-str_subset(test.features,pattern = "mm10")
test.features <- setdiff(test.features,mouse.genes)
test.features <- setdiff(test.features,featuresUCOERG)


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
  

#MouseID
DimPlot(seurat.list.act.tr.untr.MNN,group.by = "MouseID",cols = mouse.cols,pt.size = 0.5)

DimPlot(seurat.list.act.tr.untr.MNN,group.by = "MouseID",cols = mouse.cols,pt.size = 0.5)+NoLegend()

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
  geom_point()+scale_x_log10()+scale_y_log10()+labs(title = "HIV+ vs. HIV- UMIs")+
  scale_color_manual(values = HIV.cols,breaks = c("HIV+ high","HIV+ low","HIV-"))+
  theme_bw()

#Finishing pdf
dev.off()

#END OF SCRIPT####
