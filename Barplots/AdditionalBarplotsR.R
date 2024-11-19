##Skript for creating the supplementary barplots seen in figure XX and subsequent 
#This skript is applied to the final list of seurat.objects

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


#Creating table using all meta.data from the list of seurat.objects


MetadataTable<-tibble()

library(data.table)

for(i in 1:length(seurat.list)){
  metadata <- seurat.list[[i]]@meta.data
  
  # Convert metadata to data.table
  metadata_dt <- as.data.table(metadata)
  
  MetadataTable <- rbind(MetadataTable, metadata_dt)
}


#MetadataTable<-MetadataTable[MetadataTable$MouseID!="A0",]

MetadataTable2<-MetadataTable%>%group_by(MouseID)%>%dplyr::count(MouseID)

MetadataTable3<-MetadataTable%>%group_by(Count)%>%dplyr::count(Count)

MetadataTable3<-MetadataTable%>%group_by(MouseID)%>%dplyr::count(Count)

MetadataTableHIV<-MetadataTable[MetadataTable$status!="HIV-",]

MetadataTableHIV2<-MetadataTableHIV%>%group_by(MouseID)%>%dplyr::count(status)
AddHIVTable<-data.frame(c("T2","T3","U1","U2"),c("HIV+ high","HIV+ low","HIV+ high","HIV+ low","HIV+ high","HIV+ low","HIV+ high","HIV+ low"),0)#Adding values for these mice manually
AddHIVTable<-`colnames<-`(AddHIVTable,c("MouseID","status","n"))
MetadataTableHIV2<-rbind(MetadataTableHIV2,AddHIVTable)

#Fig.4E
#Creating Plots in a combined pdf format
pdf("AdditionalBarplots.pdf")

ggplot(MetadataTableHIV2,aes(x=MouseID, fill=status)) + geom_bar(position = "fill")+theme_bw()+
  theme(axis.text.x = element_text(hjust=1,angle=45),legend.key.size = unit(0.1, 'cm'))+
  scale_fill_manual(values = HIV.cols)+
  scale_y_continuous(labels = scales::percent)+
  labs(title = "# of HIV+ cells by mouse", x = "MouseID", y = "# of cells")+NoLegend()

ggplot(MetadataTableHIV2,aes(x=MouseID, fill=status,y=n)) + geom_bar(stat = "identity")+theme_bw()+
  theme(axis.text.x = element_text(hjust=1,angle=45),legend.key.size = unit(0.1, 'cm'))+
  scale_fill_manual(values = HIV.cols)+
  geom_label(aes(label=after_stat(MetadataTableHIV2$n)),position = "stack",size=3,color="black")+
  labs(title = "# of HIV+ cells by mouse", x = "MouseID", y = "# of cells")+NoLegend()

ggplot(MetadataTableHIV3,aes(x=MouseID, fill=MouseID,y=n)) + geom_bar(stat = "identity")+theme_bw()+
  theme(axis.text.x = element_text(hjust=1,angle=45),legend.key.size = unit(0.1, 'cm'))+
  scale_fill_manual(values = mouse.cols)+
  geom_label(aes(label=after_stat(MetadataTableHIV3$n)),position = "stack",size=3,color="black")+
  labs(title = "# of HIV+ cells by mouse", x = "MouseID", y = "# of cells")+NoLegend()

#Finishing the pdf
dev.off()
#END OF SCRIPT####



