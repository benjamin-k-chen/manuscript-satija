
set.seed(10403)
library(rstatix)
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggpubr)

#This is just for plotting the calculations were performed in excel
#Raw data is available in the supplementary data excel file

#Fig-4E

#Colors
HIV.cols <- c('HIV-'='lightgrey','HIV+ low'="#ffd8b1",'HIV+ high'= "purple",'HIV+ very high'="purple")
fluor.cols <- c('Unmarked'='lightgrey','Mixed'='orange','dsRed'='#f5eaea','GFP'='darkgreen')

thresh.cols <- c('N'='lightgrey','Y'='#EF3E3E')


#Figure 4H - FastMNN with all mice (A1-3,T1-3,U1-2)#######

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


####

AcuteMedian<-(sum(Acute$HIVpos)/sum(Acute$Total))*100

TreatedMedian<-(sum(Treated$GFP)/sum(Treated$Total))*100

####
#Created plot


Acute$Threshold<-case_when(Acute$Total>=200 ~ "Y",.default = "N" )

Treated$Threshold<-case_when(Treated$Total>=200 ~ "Y",.default = "N" )

P1<-ggplot(Acute,aes(x=CellType, y=`Log2FC Acute`, fill=Threshold)) + geom_col()+theme_bw()+
  scale_x_discrete(limits=rev(cell.orders))+
  scale_y_continuous(limits = c(-2.6,2.6))+
  coord_flip()+
  scale_fill_manual(values = thresh.cols)+
  NoLegend()

P2<-ggplot(Treated,aes(x=CellType, y=`Log2FC Treated`, fill=Threshold)) + geom_col()+theme_bw()+
  scale_x_discrete(limits=rev(cell.orders))+
  scale_y_continuous(limits = c(-2.6,2.6))+
  coord_flip()+
  scale_fill_manual(values = thresh.cols)+
  NoLegend()

ggarrange(P1,P2,ncol = 2,nrow = 1)

###

#Some kind of coefficient??
#LogFC Acute/Log2FC Treated

Acute$Log2FCratio <- Treated$`Log2FC Treated`-Acute$`Log2FC Acute`









