#FastMNN skript 1st try
#Setting the seed to allow reproducibility in terms of UMAPs
set.seed(10403)

library(SeuratWrappers)
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyverse)


#load the seurat list of the merges from the pipeline
seurat.list.acute<-readRDS("/Users/gerritschmidt/Desktop/seurat.list.acuteMNN.rds") #SCT transformed and merged on Minerva (simply because these are big)

#seurat.list.act.tr<-readRDS("/Users/gerritschmidt/Documents/scRNA docs/Results-by-reference/MergedByMouse/Hope/20230314-FastFNN+FullPipeline/seurat.list.act.tr.rds") #SCT transformed and merged on Minerva (simply because these are big)

#seurat.list.act.tr.uninf<-readRDS("/Users/gerritschmidt/Desktop//seurat.list.act.tr.uninf.rds") #SCT transformed and merged on Minerva (simply because these are big)

seurat.list.act.tr.uninf.MNN<-readRDS("/Users/gerritschmidt/Desktop/seurat.list.act.tr.untr.MNN.rds") #SCT transformed and merged on Minerva (simply because these are big)



featuresUCOERG<- c("dsRed", "EGFP", "WPRE" )
featuresHIV<- c("gag-pol","pol","pol-vif-vpr-tat-rev","vpu-env","env","mirfp670nano","p2a-cre-ires-nef" )
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

datasetsMouseID<-c("27M1"="T3",
              "297_count" = "A1",
              "31M1" = "A0",
              "368_count" = "T2",
              "439_count" = "T1",
              "440_count" = "T1",
              "609_count" = "U1",
              "610_count" = "U2",
              "611_count" = "U1",
              "612_count" = "U2",
              "618_count" = "A2",
              "619_count" = "A2",
              "620_count" = "A3",
              "621_count" = "A3",
              "622_count" = "A3")


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


#Add MouseID in MetaData

datasetsMetaData <- FetchData(seurat.list.act.tr.uninf.MNN,vars = c("Count","Fluorescence"),slot = "counts")
datasetsMetaData$cell.id <- rownames(datasetsMetaData)

datasetsMouseIDtbl<-data.frame(datasetsMouseID)
datasetsMouseIDtbl$Count<-rownames(datasetsMouseIDtbl)

seurat.info<-inner_join(datasetsMetaData,datasetsMouseIDtbl,by="Count")
rownames(seurat.info)<-seurat.info$cell.id

seurat.list.act.tr.uninf.MNN <- AddMetaData(seurat.list.act.tr.uninf.MNN, metadata = seurat.info$datasetsMouseID, col.name = "MouseID")


#Rename predicted.celltypes "CD8" => "CD8 Naive" "CD8 Proliferating" "CD8 TEM"
# "Others" -> "CD14 Mono" , "cDC1" , "cDC2", "HSPC", "Platelet"

seurat.list.acute.MNN.copy<-seurat.list.acute

seurat.list.acute$predicted.celltype.l2old<-seurat.list.acute$predicted.celltype.l2
seurat.list.acute$predicted.celltype.l2 <- str_replace_all(seurat.list.acute$predicted.celltype.l2,c("CD8 Naive|CD8 Proliferating|CD8 TEM"),"CD8")
seurat.list.acute$predicted.celltype.l2 <- str_replace_all(seurat.list.acute$predicted.celltype.l2,c("CD14 Mono|cDC1|cDC2|HSPC|Platelet"),"Others")

######
##Acute-TR-Uninf

#test.features2 <- VariableFeatures(seurat.list.act.tr.uninf)
#test.features2 <- setdiff(test.features2,featuresHIV)
#mouse.genes<-str_subset(test.features2,pattern = "mm10")
#test.features2 <- setdiff(test.features2,mouse.genes)
#test.features2 <- setdiff(test.features2,featuresUCOERG)


#seurat.list.act.tr.uninf <- RunFastMNN(object.list = SplitObject(seurat.list.act.tr.uninf,split.by = "Count"), features=test.features2)

seurat.list.act.tr.uninf.MNN <- RunUMAP(seurat.list.act.tr.uninf.MNN,dims = 1:30,reduction = "mnn")


#
cell.order.abundance<-seurat.list.act.tr.uninf.MNN@meta.data%>%
   dplyr::count(predicted.celltype.l2)%>%
   mutate(percent=(n/sum(n))*100)

cell.order.abundance<-arrange(cell.order.abundance,by_group=percent)

cell.order.abundance<-cell.order.abundance$predicted.celltype.l2

##
pdf("20230424-FastMNN-AcuteMNN.pdf",width = 14,height = 10)
DimPlot(seurat.list.act.tr.uninf.MNN,group.by = "Condition",order = c("Acute","ART treated","Uninfected"),pt.size = 0.5)+
  scale_color_manual(values = cond.cols,breaks = c("Acute","ART treated","Uninfected"))

DimPlot(seurat.list.act.tr.uninf.MNN,group.by = "Condition",order = c("Acute","ART treated","Uninfected"),pt.size = 0.5)+
  scale_color_manual(values = cond.cols,breaks = c("Acute","ART treated","Uninfected"))+NoLegend()
  
DimPlot(seurat.list.act.tr.uninf.MNN,group.by = "MouseID",cols = mouse.cols,pt.size = 0.5)

DimPlot(seurat.list.act.tr.uninf.MNN,group.by = "MouseID",cols = mouse.cols,pt.size = 0.5)+NoLegend()


DimPlot(seurat.list.act.tr.uninf.MNN,group.by = "status",order = c("HIV+ high","HIV+ low","HIV-"),pt.size = 0.5)+
  scale_color_manual(values = HIV.cols,breaks = c("HIV+ high","HIV+ low","HIV-"))

DimPlot(seurat.list.act.tr.uninf.MNN,group.by = "status",order = c("HIV+ high","HIV+ low","HIV-"),pt.size = 0.5)+
  scale_color_manual(values = HIV.cols,breaks = c("HIV+ high","HIV+ low","HIV-"))+NoLegend()

DimPlot(seurat.list.act.tr.uninf.MNN,group.by = "predicted.celltype.l2",order = cell.order.abundance,pt.size = 0.5)+
  scale_color_manual(values = cell.cols,breaks = cell.order.rev)

DimPlot(seurat.list.act.tr.uninf.MNN,group.by = "predicted.celltype.l2",order = cell.order.abundance,pt.size = 0.5)+
  scale_color_manual(values = cell.cols,breaks = cell.order.rev)

DimPlot(seurat.list.act.tr.uninf.MNN,group.by = "Fluorescence",pt.size = 0.5,cols = fluor.cols,order = c("GFP","dsRed","Mixed","Unmarked"))
DimPlot(seurat.list.act.tr.uninf.MNN,group.by = "Fluorescence",pt.size = 0.5,cols = fluor.cols,order = c("GFP","dsRed","Mixed","Unmarked"))+
  NoLegend()

DimPlot(seurat.list.act.tr.uninf.MNN,group.by = "Doublet.Singlet",pt.size = 0.5,order = c("Doublet","Singlet"),cols = c("Doublet"="darkblue","Singlet"="lightgrey"))
DimPlot(seurat.list.act.tr.uninf.MNN,group.by = "Doublet.Singlet",pt.size = 0.5,order = c("Doublet","Singlet"),cols = c("Doublet"="darkblue","Singlet"="lightgrey"))+
  NoLegend()


##Combined UMI plot

seurat.list.act.tr.uninf.MNN.AandT<-subset(seurat.list.act.tr.uninf.MNN,subset = MouseID %in% c("A0","A1","A2","A3","T1","T2","T3"))

gene_HIV2 <- ifelse(str_detect(rownames(seurat.list.act.tr.uninf.MNN.AandT@assays$RNA),"gag-pol|^pol|pol-vif-vpr-tat-rev|vpu-env|env|p2a-cre-ires-nef"),"HIV+", "HIV-")
hivpos_inds2 <- gene_HIV2 == "HIV+"
hivneg_inds2 <- gene_HIV2 == "HIV-"
cell_hivstat2 <- tibble(n_hivpos_umi = Matrix::colSums(seurat.list.act.tr.uninf.MNN.AandT@assays$RNA[hivpos_inds2,]),
                        n_hivneg_umi = Matrix::colSums(seurat.list.act.tr.uninf.MNN.AandT@assays$RNA[hivneg_inds2,]),
                        tot_umi = Matrix::colSums(seurat.list.act.tr.uninf.MNN.AandT@assays$RNA),
                        Status = seurat.list.act.tr.uninf.MNN.AandT$status,
                        Count = seurat.list.act.tr.uninf.MNN.AandT$Count,
                        Condition = seurat.list.act.tr.uninf.MNN.AandT$Condition,
                        MouseID = seurat.list.act.tr.uninf.MNN.AandT$MouseID,
                        Mean_hivpos_umi = base::mean(n_hivpos_umi))
cell_hivstat2<-cell_hivstat2%>%
  arrange(Status)
ggplot(cell_hivstat2, aes(n_hivpos_umi, n_hivneg_umi, color = Status,shape=MouseID)) +
  geom_point()+scale_x_log10()+scale_y_log10()+labs(title = "HIV+ vs. HIV- UMIs")+
  scale_color_manual(values = HIV.cols,breaks = c("HIV+ high","HIV+ low","HIV-"))+
  theme_bw()

dev.off()



ggplot(seurat.list.act.tr.uninf.MNN@meta.data,aes(x=seurat_clusters, fill=Predicted.Celltype.L2)) + geom_bar(position="fill")+
  scale_y_continuous(labels = scales::percent)+theme_bw()+labs(title = paste0(dplyr::first(seurat.list.act.tr.uninf.MNN$Count)),
                                                               subtitle = "After cluster subsetting")+
  theme(axis.text.x = element_text(hjust=1,angle=45),legend.key.size = unit(0.1, 'cm'))+
  scale_fill_manual(values = cell.cols,breaks = cell.order.rev)


cell_table <- seurat.list.act.tr.uninf.MNN@meta.data%>%
  group_by(seurat_clusters)%>%
  dplyr::count(predicted.celltype.l2)%>%
  mutate(percent=(n/sum(n))*100)

cell_table2<-cell_table[order(match(cell_table$predicted.celltype.l2, cell.order.rev)),]

ggplot(cell_table2,aes(x=seurat_clusters,fill=predicted.celltype.l2,y=percent)) + geom_col(position = position_stack(reverse = FALSE))+
  theme_bw()+labs(title = paste0(dplyr::first(seurat.list.act.tr.uninf.MNN$Count)),
                                                               subtitle = "After cluster subsetting")+
  theme(axis.text.x = element_text(hjust=1,angle=45),legend.key.size = unit(0.1, 'cm'))+
  scale_fill_manual(values = cell.cols,breaks = cell.order.rev)


limits = c(0,0.20)





#####
#AcuteMerge
#seurat.list.acute<-merge(`297_count`,c(`618_count`,`620_count`))

#seurat.list.acute<-SCTransform(seurat.list.acute)

test.features <- VariableFeatures(seurat)
test.features <- setdiff(test.features,featuresHIV)
mouse.genes<-str_subset(test.features,pattern = "mm10")
test.features <- setdiff(test.features,mouse.genes)
test.features <- setdiff(test.features,featuresUCOERG)

##https://github.com/satijalab/seurat/issues/5329 (The error that RunFastMNN will throw)
# this is how they say to fix it
trace('RunFastMNN', edit = T)
# change line 25-28
# objects.sce <- lapply(X = object.list, FUN = function(x,f) { if (DefaultAssay(x) == 'SCT') { x = subset(x = x, features = f) indx <- match(rownames(x@assays$SCT@counts),rownames(x@assays$SCT@scale.data)) x@assays$SCT@scale.data <- x@assays$SCT@scale.data[indx,] }else{ x = subset(x = x, features = f) } return(as.SingleCellExperiment(x)) }, f = features)
#
#

seurat.list.acute <- RunFastMNN(object.list = SplitObject(seurat.list.acute,split.by = "Count"), features=test.features)

seurat.list.acute <- RunUMAP(seurat.list.acute,dims = 1:30,reduction = "mnn")

acutep1<-DimPlot(seurat.list.acute,group.by = "Treatment")
acutep2<-DimPlot(seurat.list.acute,group.by = "SampleType")

#seurat.list.acute <- FindNeighbors(seurat.list.acute,dims = 1:unique(seurat.info$datasetsUMAP))
#seurat.list.acute <- FindClusters(seurat.list.merged.by.mouse.copy.subset[[i]],resolution = resolution.range) 


test.features <- VariableFeatures(seurat.list.act.tr)
test.features <- setdiff(test.features,featuresHIV)
mouse.genes<-str_subset(test.features,pattern = "mm10")
test.features <- setdiff(test.features,mouse.genes)
test.features <- setdiff(test.features,featuresUCOERG)


seurat.list.act.tr <- RunFastMNN(object.list = SplitObject(seurat.list.act.tr,split.by = "Count"), features=test.features)

seurat.list.act.tr <- RunUMAP(seurat.list.act.tr,dims = 1:30,reduction = "mnn")

acutetrp3<-DimPlot(seurat.list.act.tr,group.by = "Treatment")
acutetrp4<-DimPlot(seurat.list.act.tr,group.by = "SampleType")

substr(Sys.time(),1,10)
write_csv(table1stUmap618,paste0("Table618_1stUMAP-",substr(Sys.time(),1,10),".csv"))

