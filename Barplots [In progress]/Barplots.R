##Skript for creating the distribution barplots seen in figure XX and subsequent individual supplementary figures
#This skript is applied to the final list of seurat objects

library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggpubr)


#We first load our final seurat object list
seurat.list<-readRDS("seurat.rds")

#Colors and legend order#####
#We then define colors, genes and legend order which we will need for finetuning the later plots
#These are the same as in the main script as well in the FastMNN script

dsRed.cols <- c('LentiRG'='red','No expression'='lightgrey')
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

#Creating output pdf file####
pdf("Barplots.pdf")

#Loop for plots#####
for (i in 1:length(seurat.list)){
  
  #Barplots - predicted.celltype.l1 and l2
  
  #Creating a celltype.l2 table
  cell_table <- seurat.list[[i]]@meta.data%>%
    group_by(seurat_clusters)%>%
    dplyr::count(predicted.celltype.l2)%>%
    mutate(percent=(n/sum(n)))
  
  cell_table2<-cell_table[order(match(cell_table$predicted.celltype.l2, cell.order.rev)),]
  
  cell_table2$predicted.celltype.l2<-factor(cell_table2$predicted.celltype.l2,levels = unique(cell_table2$predicted.celltype.l2))
  
  a1<-ggplot(cell_table2,aes(x=seurat_clusters,y=cell_table2$percent,fill=factor(cell_table2$predicted.celltype.l2,levels = cell.order.rev)))+
    geom_bar(stat = "identity",position = "stack")+
    scale_y_continuous(expand=c(0,0),labels = scales::percent)+
    theme_bw()+labs(fill="predicted.celltype.l2",y="% of cells",x="Cluster",caption = paste0(dplyr::first(seurat.list[[i]]$Count),"-final barplots"))+
    theme(axis.text.x = element_text(hjust=1,angle=45),legend.key.size = unit(0.1, 'cm'))+
    scale_fill_manual(values = cell.cols,breaks = cell.order.rev)
  
  #cerating the .l1 barplots (-> celltypes are ordered alphabetically)
  a2<-ggplot(seurat.list[[i]]@meta.data,aes(x=seurat_clusters, fill=predicted.celltype.l1))+
    geom_bar(position="fill")+
    scale_y_continuous(labels = scales::percent,expand=c(0,0))+
    theme_bw()+
    labs(fill="predicted.celltype.l1",y="% of cells",x="Cluster",caption = paste0(dplyr::first(seurat.list[[i]]$Count),"-final barplots"))+
    theme(axis.text.x = element_text(hjust=1,angle=45),legend.key.size = unit(0.1, 'cm'))

  show(ggarrange(a1,a2,ncol=1,nrow = 2))
  
  #Barplots -HIV status
  
  h1<-ggplot(seurat.list[[i]]@meta.data,aes(x=seurat_clusters, fill=status)) + geom_bar(position="fill")+
    scale_y_continuous(labels = scales::percent,expand=c(0,0))+
    theme_bw()+
    labs(fill="HIV status",y="% of cells",x="Cluster",caption = paste0(dplyr::first(seurat.list[[i]]$Count),"-final barplots"))+
    theme(axis.text.x = element_text(hjust=1,angle=45),legend.key.size = unit(0.1, 'cm'))+
    scale_fill_manual(values = HIV.cols,breaks = c("HIV+ high","HIV+ low","HIV-"))+
    coord_cartesian(ylim = c(0,0.25))
  
  h2<-ggplot(seurat.list[[i]]@meta.data,aes(x=predicted.celltype.l2, fill=status)) + geom_bar(position="fill")+
    scale_y_continuous(labels = scales::percent,expand=c(0,0))+
    theme_bw()+
    labs(fill="HIV status",y="% of cells",x="predicted.celltypes.l2",caption = paste0(dplyr::first(seurat.list[[i]]$Count),"-final barplots"))+
    theme(axis.text.x = element_text(hjust=1,angle=45),legend.key.size = unit(0.1, 'cm'))+
    scale_fill_manual(values = HIV.cols,breaks = c("HIV+ high","HIV+ low","HIV-"))+
    coord_cartesian(ylim = c(0,0.25))
  show(ggarrange(h1,h2,ncol=1,nrow = 2))

  #Barplots - LentiRG expression
  
  t1<-ggplot(seurat.list[[i]]@meta.data,aes(x=seurat_clusters, fill=LentiRG.expression))+geom_bar(position="fill")+
    scale_y_continuous(labels = scales::percent,expand=c(0,0))+theme_bw()+
    labs(fill="LentiRG status",y="% of cells",x="Cluster",caption = paste0(dplyr::first(seurat.list[[i]]$Count),"-final barplots"))+
    theme(axis.text.x = element_text(hjust=1,angle=45),legend.key.size = unit(0.1, 'cm'))+
    scale_fill_manual(values = dsRed.cols)
  
  t2<-ggplot(seurat.list[[i]]@meta.data,aes(x=predicted.celltype.l2, fill=LentiRG.expression))+geom_bar(position="fill")+
    scale_y_continuous(labels = scales::percent,expand=c(0,0))+theme_bw()+
    labs(fill="LentiRG status",y="% of cells",x="predicted.celltypes.l2",caption = paste0(dplyr::first(seurat.list[[i]]$Count),"-final barplots"))+
    theme(axis.text.x = element_text(hjust=1,angle=45),legend.key.size = unit(0.1, 'cm'))+
    scale_fill_manual(values = dsRed.cols)
  
  show(ggarrange(t1,t2,ncol=1,nrow = 2))
  
  #Barplots - Fluorescence
  o1<-ggplot(seurat.list[[i]]@meta.data,aes(x=seurat_clusters, fill=Fluorescence))+geom_bar(position="fill")+
    scale_y_continuous(labels = scales::percent,expand=c(0,0))+theme_bw()+
    labs(fill="Fluorescence",y="% of cells",x="Cluster",caption = paste0(dplyr::first(seurat.list[[i]]$Count),"-final barplots"))+
    theme(axis.text.x = element_text(hjust=1,angle=45),legend.key.size = unit(0.1, 'cm'))+
    scale_fill_manual(values = fluor.cols)
  
  o2<-ggplot(seurat.list[[i]]@meta.data,aes(x=predicted.celltype.l2, fill=Fluorescence)) + geom_bar(position="fill")+
    scale_y_continuous(labels = scales::percent,expand=c(0,0))+theme_bw()+
    labs(fill="Fluorescence",y="% of cells",x="predicted.celltype.l2",caption = paste0(dplyr::first(seurat.list[[i]]$Count),"-final barplots"))+
    theme(axis.text.x = element_text(hjust=1,angle=45),legend.key.size = unit(0.1, 'cm'))+
    scale_fill_manual(values = fluor.cols)
  
  show(ggarrange(o1,o2,ncol=1,nrow = 2))
}

#Finishing .pdf file####  
dev.off()

#END OF SCRIPT####






