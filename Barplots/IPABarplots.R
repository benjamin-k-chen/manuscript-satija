#IPA_Barplots
##Skript for creating the IPA barplots seen in figure XX and subsequent individual supplementary figures
#This skript uses the outpot .csv pathway list from IPA

library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggpubr)

#Fig. 5B-C-E-F
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


#Barplots using Enrichr Results
#This is just the plotting function based on the obtained Enrichr results
#Raw data is also available in the supplementary data file

pdf("Enrichr_pvalue_barplots.pdf",width = 10)
E1<-ggplot()+geom_bar(T1T2_Enrichr,stat = "identity",mapping=aes(x=-log(`P-value`),y=reorder(Pathways,-log(`P-value`)),fill=-log(`P-value`)))+
  scale_fill_gradient2(high = "maroon",limits=c(0,30))+
  theme_pubr()+labs(title = "Enrichr - Analysis - T1T2")+xlab("-log(p-value)")+ylab("")

E2<-ggplot()+geom_bar(A3_Enrichr,stat = "identity",mapping=aes(x=-log(`P-value`),y=reorder(Pathways,-log(`P-value`)),fill=-log(`P-value`)))+
  scale_fill_gradient2(high = "maroon",limits=c(0,30))+
  theme_pubr()+labs(title = "Enrichr - Analysis - A3 Acute")+xlab("-log(p-value)")+ylab("")
ggarrange(E1,E2,align = "hv",ncol = 1)
dev.off()




##Fig. 5G-H
pdf("Alternative_IPA_Comparison_analysis.pdf",width = 20)

PathwaysPValue<-Satija_Acute_celltypes_vs_Treated_pvalue$`Canonical Pathways`[1:15]
P1<-ggplot(data = Satija_Acute_celltypes_vs_Treated_pvalue, aes(y = `Canonical Pathways`, x = celltype)) +
  geom_tile(aes(fill = p_value)) +
  scale_x_discrete(position = "top")+
  ylim(rev(PathwaysPValue))+
  scale_fill_gradient2(high = "purple",low = "lightgrey",mid = "white",na.value = "#e5e5e5")+
  theme_minimal()


PathwaysZScore<-Satija_Acute_celltypes_vs_Treated_zscore$`Canonical Pathways`[1:15]
P2<-ggplot(data = Satija_Acute_celltypes_vs_Treated_zscore, aes(y = `Canonical Pathways`, x = celltype)) +
  geom_tile(aes(fill = z_score)) +
  scale_x_discrete(position = "top")+
  ylim(rev(PathwaysZScore))+
  scale_fill_gradient2(high = "orange",low = "blue",mid = "white",na.value = "#e5e5e5")+
  theme_minimal()

ggarrange(P1,P2,align = "hv")

dev.off()


##Fig. 6A-B
pdf("Alternative_IPA_Comparison_analysis.pdf",width = 20)

PathwaysPValue<-Collora_vs_Satija_pvalue$`Canonical Pathways`[1:15]
P1<-ggplot(data = Collora_vs_Satija_pvalue, aes(y = `Canonical Pathways`, x = Dataset)) +
  geom_tile(aes(fill = p_value)) +
  scale_x_discrete(position = "top")+
  ylim(rev(PathwaysPValue))+
  scale_fill_gradient2(high = "purple",low = "lightgrey",mid = "white",na.value = "#e5e5e5")+
  theme_minimal()


PathwaysZScore<-Collora_vs_Satija_zscore$`Canonical Pathways`[1:15]
P2<-ggplot(data = Collora_vs_Satija_zscore, aes(y = `Canonical Pathways`, x = Dataset)) +
  geom_tile(aes(fill = z_score)) +
  scale_x_discrete(position = "top")+
  ylim(rev(PathwaysZScore))+
  scale_fill_gradient2(high = "orange",low = "blue",mid = "white",na.value = "#e5e5e5")+
  theme_minimal()

ggarrange(P1,P2,align = "hv")

dev.off()


##Fig. 6C - Fig.S7A
pdf("Alternative_IPA_Comparison_analysis.pdf",width = 20)

PathwaysPValue<-Acute_vs_Treated_pvalue$`Canonical Pathways`[1:15]
P1<-ggplot(data = Acute_vs_Treated_pvalue, aes(y = `Canonical Pathways`, x = Mouse)) +
  geom_tile(aes(fill = p_value)) +
  scale_x_discrete(position = "top")+
  ylim(rev(PathwaysPValue))+
  scale_fill_gradient2(high = "purple",low = "lightgrey",mid = "white",na.value = "#e5e5e5")+
  theme_minimal()


PathwaysZScore<-Acute_vs_Treated_zscore$`Canonical Pathways`[1:15]
P2<-ggplot(data = Acute_vs_Treated_zscore, aes(y = `Canonical Pathways`, x = Mouse)) +
  geom_tile(aes(fill = z_score)) +
  scale_x_discrete(position = "top")+
  ylim(rev(PathwaysZScore))+
  scale_fill_gradient2(high = "orange",low = "blue",mid = "white",na.value = "#e5e5e5")+
  theme_minimal()

ggarrange(P1,P2,align = "hv")

dev.off()