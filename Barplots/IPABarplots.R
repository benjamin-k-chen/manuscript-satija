#IPA_Barplots
##Skript for creating the IPA barplots seen in figure XX and subsequent individual supplementary figures
#This skript uses the outpot .csv pathway list from IPA

library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggpubr)

pdf("Alternative_IPA_pathway_analysis.pdf",width = 20)
A1_CannonicalPathwaysShort<-A1_CannonicalPathways[c(1:5),]
p1<-ggplot()+geom_bar(A1_CannonicalPathwaysShort,stat = "identity",mapping=aes(x=`-log(p-value)`,y=reorder(`Ingenuity Canonical Pathways`,`-log(p-value)`),fill=`z-score`))+
  scale_fill_gradient2(high = "orange",low = "blue",mid = "lightgrey",na.value = "#e5e5e5")+
  theme_pubr()+labs(title = "IPA - Analysis - A1")+xlab("-log(p-value)")+ylab("")

A2_CannonicalPathwaysShort<-A2_CannonicalPathways[c(1:5),]
p2<-ggplot()+geom_bar(A2_CannonicalPathwaysShort,stat = "identity",mapping=aes(x=`-log(p-value)`,y=reorder(`Ingenuity Canonical Pathways`,`-log(p-value)`),fill=`z-score`))+
  scale_fill_gradient2(high = "orange",low = "blue",mid = "lightgrey",na.value = "#e5e5e5")+
  theme_pubr()+labs(title = "IPA - Analysis - A2")+xlab("-log(p-value)")+ylab("")

A3_CannonicalPathwaysShort<-A3_CannonicalPathways[c(1:5),]
p3<-ggplot()+geom_bar(A3_CannonicalPathwaysShort,stat = "identity",mapping=aes(x=`-log(p-value)`,y=reorder(`Ingenuity Canonical Pathways`,`-log(p-value)`),fill=`z-score`))+
  scale_fill_gradient2(high = "orange",low = "blue",mid = "lightgrey",na.value = "#e5e5e5")+
  theme_pubr()+labs(title = "IPA - Analysis - A3")+xlab("-log(p-value)")+ylab("")

T1T2_CannonicalPathwaysShort<-T1T2_CannonicalPathways[c(1:5),]
p4<-ggplot()+geom_bar(T1T2_CannonicalPathwaysShort,stat = "identity",mapping=aes(x=`-log(p-value)`,y=reorder(`Ingenuity Canonical Pathways`,`-log(p-value)`),fill=`z-score`))+
  scale_fill_gradient2(high = "orange",low = "blue",mid = "lightgrey",na.value = "#e5e5e5")+
  theme_pubr()+labs(title = "IPA - Analysis - T1T2")+xlab("-log(p-value)")+ylab("")
ggarrange(p1,p2,p3,p4,align = "hv")
dev.off()