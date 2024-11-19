##Lets add Fisher exact test to selected barplots

set.seed(10403)
library(rstatix)
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggpubr)

#pairwise_fisher_test():Performs Fisher's exact test for testing the null of independence of rows and columns in a contingency table.
#Wrappers around the R base function fisher.test() but have the advantage of performing pairwise and row-wise fisher tests, 
#the post-hoc tests following a significant chi-square test of homogeneity for 2xc and rx2 contingency tables.


#Utilities

#Colors
HIV.cols <- c('HIV-'='lightgrey','HIV+ low'="#ffd8b1",'HIV+ high'= "purple",'HIV+ very high'="purple")
fluor.cols <- c('Unmarked'='lightgrey','Mixed'='orange','dsRed'='#f5eaea','GFP'='darkgreen')

#Figure 3H - FastMNN with all mice (A1-3,T1-3,U1-2)#######

cell.orders<-c("Proliferating CD4",
               "Th1 Proliferating CD4",
               "Cytolytic Proliferating CD4",
               "Naïve CD4-1",
               "Naïve CD4-2",
               "Cytolytic CD4 TCM",
               "Tfh TCM",
               "Th1/Th17 TEM",
               "Th22/Th9 TEM",
               "Tfh/Th22",
               "Th2",
               "Tregs",
               "IL10RA+ Tregs",
               "MT Cells")

cell.orders.clusters<-c("1","5","11","8","9","10","2","0","13","6","12","3","4","7")

#Let's extract the numbers of cells HIV-/HIV+ using cluster designations

total_vector<-seurat.list.act.mtr.untr.MNN.clusters@meta.data%>%
  group_by(seurat_clusters)%>%
  dplyr::count(status=="HIV-")

HIVpos.table<-total_vector[total_vector$`status == "HIV-"`=="FALSE",]
HIVneg.table<-total_vector[total_vector$`status == "HIV-"`=="TRUE",]

total_vector1<-seurat.list.act.mtr.untr.MNN.clusters@meta.data%>%
  group_by(seurat_clusters)%>%
  dplyr::count(seurat_clusters)

TableA <- data.frame(
  HIVpos=c(total_vector1$n-HIVneg.table$n),#Sucesses
  HIVneg=HIVneg.table$n#Failures
)

#Unfortunately I'm using a version of the seurat object which doesn't have the final cluster desigantions so I'll add them manually
TableA<-`rownames<-`(TableA,c("Th1/Th17 TEM",
                              "Proliferating CD4",
                              "Tfh TCM",
                              "Tregs",
                              "IL10RA+ Tregs",
                              "Th1 Proliferating CD4",
                              "Tfh/Th22",
                              "MT Cells",
                              "Naïve CD4-1",
                              "Naïve CD4-2",
                              "Cytolytic CD4 TCM",
                              "Cytolytic Proliferating CD4",
                              "Th2",
                              "Th22/Th9 TEM"
))

##Lets do the Fisher's Exact Test

#Fishers Exact Test
TableA_fisher<-fisher_test(TableA,simulate.p.value = TRUE)

#Fishers Exact Test - Post hoc pairwise
pairwisetest_Fig3H<-pairwise_fisher_test(x=TableA,p.adjust.method="bonferroni")

#Comparing each cluster to the overall average

TableA_results<-tibble()
for(i in 1:nrow(TableA)) {
  TableA_Total<-data_frame(HIVpos=sum(TableA$HIVpos),HIVneg=sum(TableA$HIVneg))
  if (TableA$HIVpos[i]>1){
    print("Doing Fishers exact test")
    TableA1<-TableA[i,]
    TableA1<-rbind(TableA1,TableA_Total)
    res<-pairwise_fisher_test(TableA1,p.adjust.method = "bonferroni")
    print(res)
    TableA_results<-rbind(TableA_results,res)
  }
  else{
    print("No HIV+ -> skip")
  }
}

####Create two barplot versions
a<-ggplot(seurat.list.act.mtr.untr.MNN.clusters@meta.data,aes(x=seurat_clusters, fill=status))+
  geom_bar(position="fill")+
  geom_hline(yintercept = TableA_Total$HIVpos/(TableA_Total$HIVneg+TableA_Total$HIVpos), color="darkgrey",linetype=2)+
  scale_y_continuous(labels = scales::percent)+
  scale_x_discrete(limits=cell.orders.clusters,labels=cell.orders)+
  theme_bw()+
  labs(fill="HIV status",y="% of cells",x="")+
  theme(axis.text.x = element_text(hjust=1,angle=45),legend.key.size = unit(0.1, 'cm'))+
  coord_cartesian(ylim = c(0,0.05))+
  scale_fill_manual(values = HIV.cols,breaks = c("HIV+ high","HIV+ low","HIV-"))

b<-ggplot(seurat.list.act.mtr.untr.MNN.clusters@meta.data,aes(x=seurat_clusters, fill=status))+
  geom_bar(stat ="count",position="stack")+
  scale_x_discrete(limits=cell.orders.clusters,labels=cell.orders)+
  theme_bw()+
  labs(fill="HIV status",y="# of cells",x="")+
  theme(axis.text.x = element_text(hjust=1,angle=45),legend.key.size = unit(0.1, 'cm'))+
  scale_fill_manual(values = HIV.cols,breaks = c("HIV+ high","HIV+ low","HIV-"))

####Create heatmaps to visualize pairwise results

c<-ggplot(pairwisetest_Fig3H,
       aes(group1, group2))+ # x and y axes => Var1 and Var2
  geom_tile(aes(fill = p.adj))+ # background colours are mapped according to the value column
  geom_text(aes(label = p.adj),size=2)+ # write the values
  scale_fill_gradient2(low = "darkred",
                       mid = "white",
                       high = "white",
                       midpoint=0.05,
                       limits=c(min(pairwisetest_Fig3H$p.adj),max(pairwisetest_Fig3H$p.adj)),
                       trans="log") + # determine the colour
  scale_x_discrete(limits=unique(pairwisetest_Fig3H$group1),name="")+
  scale_y_discrete(limits=unique(pairwisetest_Fig3H$group2),name="")+
  NoLegend()+
  theme(panel.grid.major.x=element_blank(), 
        panel.grid.minor.x=element_blank(), 
        panel.grid.major.y=element_blank(), 
        panel.grid.minor.y=element_blank(),
        panel.background=element_rect(fill="white"),
        axis.text.x = element_text(angle=45, hjust = 1,vjust=1,size = 8),
        axis.text.y = element_text(size = 8))

d<-ggplot(pairwisetest_Fig3H,
       aes(group1, group2))+ # x and y axes => Var1 and Var2
  geom_tile(aes(fill = p.adj.signif))+ # background colours are mapped according to the value column
  geom_text(aes(label = p.adj.signif))+ # write the values
  scale_fill_manual(values=c("****" = "darkred",
                             "***" = "red",
                             "**"="orange",
                             "*"="lightyellow",
                             "ns" = "white"))+ # determine the colour
  scale_x_discrete(limits=unique(pairwisetest_Fig3H$group1),name="")+
  scale_y_discrete(limits=unique(pairwisetest_Fig3H$group2),name="")+
  NoLegend()+
  theme(panel.grid.major.x=element_blank(), 
        panel.grid.minor.x=element_blank(), 
        panel.grid.major.y=element_blank(), 
        panel.grid.minor.y=element_blank(),
        panel.background=element_rect(fill="white"),
        axis.text.x = element_text(angle=45, hjust = 1,vjust=1,size = 8),
        axis.text.y = element_text(size = 8))

pdf("Statistics_Fig3H.pdf",height = 14)
ggarrange(a,b,ncol = 1,nrow = 2)
ggarrange(c,d,ncol = 1,nrow = 2)
dev.off()


#Figure S.Fig4C - FastMNN with all mice (A1-3,T1-3,U1-2)-predicted celltype#######

#Let's extract the numbers of cells HIV-/HIV+ using cluster designations

total_vector<-seurat.list.act.mtr.untr.MNN.clusters@meta.data%>%
  group_by(predicted.celltype.l2)%>%
  dplyr::count(status=="HIV-")

HIVpos.table<-total_vector[total_vector$`status == "HIV-"`=="FALSE",]
HIVneg.table<-total_vector[total_vector$`status == "HIV-"`=="TRUE",]

total_vector1<-seurat.list.act.mtr.untr.MNN.clusters@meta.data%>%
  group_by(predicted.celltype.l2)%>%
  dplyr::count(predicted.celltype.l2)

TableB <- data.frame(
  HIVpos=c(total_vector1$n-HIVneg.table$n),#Sucesses
  HIVneg=HIVneg.table$n#Failures
)

#Unfortunately I'm using a version of the seurat object which doesn't have the final cluster desigantions so I'll add them manually
TableB<-`rownames<-`(TableB,HIVneg.table$predicted.celltype.l2)

##Lets do the Fisher's Exact Test

#Fishers Exact Test
TableB_fisher<-fisher_test(TableB,simulate.p.value = TRUE)

#Fishers Exact Test - Post hoc pairwise
pairwisetest_S.Fig4D<-pairwise_fisher_test(x=TableB,p.adjust.method="bonferroni")

#Comparing each cluster to the overall average

TableB_results<-tibble()
for(i in 1:nrow(TableB)) {
  TableB_Total<-data_frame(HIVpos=sum(TableB$HIVpos),HIVneg=sum(TableB$HIVneg))
  if (TableB$HIVpos[i]>1){
    print("Doing Fishers exact test")
    TableB1<-TableB[i,]
    TableB1<-rbind(TableB1,TableB_Total)
    res<-pairwise_fisher_test(TableB1,p.adjust.method = "bonferroni")
    print(res)
    TableB_results<-rbind(TableB_results,res)
  }
  else{
    print("No HIV+ -> skip")
  }
}

####Create two barplot versions
a<-ggplot(seurat.list.act.mtr.untr.MNN.clusters@meta.data,aes(x=predicted.celltype.l2, fill=status))+
  geom_bar(position="fill")+
  geom_hline(yintercept = TableB_Total$HIVpos/(TableB_Total$HIVneg+TableB_Total$HIVpos), color="darkgrey",linetype=2)+
  scale_y_continuous(labels = scales::percent,expand=c(0,0))+
  theme_bw()+
  labs(fill="HIV status",y="% of cells",x="")+
  theme(axis.text.x = element_text(hjust=1,angle=45),legend.key.size = unit(0.1, 'cm'))+
  scale_fill_manual(values = HIV.cols,breaks = c("HIV+ high","HIV+ low","HIV-"))+
  coord_cartesian(ylim = c(0,0.05))

b<-ggplot(seurat.list.act.mtr.untr.MNN.clusters@meta.data,aes(x=predicted.celltype.l2, fill=status))+
  geom_bar(stat ="count",position="stack")+
  theme_bw()+
  labs(fill="HIV status",y="# of cells",x="")+
  theme(axis.text.x = element_text(hjust=1,angle=45),legend.key.size = unit(0.1, 'cm'))+
  scale_fill_manual(values = HIV.cols,breaks = c("HIV+ high","HIV+ low","HIV-"))

####Create heatmaps to visualize pairwise results

c<-ggplot(pairwisetest_S.Fig4D,
          aes(group1, group2))+ # x and y axes => Var1 and Var2
  geom_tile(aes(fill = p.adj))+ # background colours are mapped according to the value column
  geom_text(aes(label = p.adj),size=2)+ # write the values
  scale_fill_gradient2(low = "darkred",
                       mid = "white",
                       high = "white",
                       midpoint=0.05,
                       limits=c(min(pairwisetest_S.Fig4D$p.adj),max(pairwisetest_S.Fig4D$p.adj)),
                       trans="log") + # determine the colour
  scale_x_discrete(limits=unique(pairwisetest_S.Fig4D$group1),name="")+
  scale_y_discrete(limits=unique(pairwisetest_S.Fig4D$group2),name="")+
  NoLegend()+
  theme(panel.grid.major.x=element_blank(), 
        panel.grid.minor.x=element_blank(), 
        panel.grid.major.y=element_blank(), 
        panel.grid.minor.y=element_blank(),
        panel.background=element_rect(fill="white"),
        axis.text.x = element_text(angle=45, hjust = 1,vjust=1,size = 8),
        axis.text.y = element_text(size = 8))

d<-ggplot(pairwisetest_S.Fig4D,
          aes(group1, group2))+ # x and y axes => Var1 and Var2
  geom_tile(aes(fill = p.adj.signif))+ # background colours are mapped according to the value column
  geom_text(aes(label = p.adj.signif))+ # write the values
  scale_fill_manual(values=c("****" = "darkred",
                             "***" = "red",
                             "**"="orange",
                             "*"="lightyellow",
                             "ns" = "white"))+ # determine the colour
  scale_x_discrete(limits=unique(pairwisetest_S.Fig4D$group1),name="")+
  scale_y_discrete(limits=unique(pairwisetest_S.Fig4D$group2),name="")+
  NoLegend()+
  theme(panel.grid.major.x=element_blank(), 
        panel.grid.minor.x=element_blank(), 
        panel.grid.major.y=element_blank(), 
        panel.grid.minor.y=element_blank(),
        panel.background=element_rect(fill="white"),
        axis.text.x = element_text(angle=45, hjust = 1,vjust=1,size = 8),
        axis.text.y = element_text(size = 8))

pdf("Statistics_SFig4D.pdf",height = 14)
ggarrange(a,b,ncol = 1,nrow = 2)
ggarrange(c,d,ncol = 1,nrow = 2)
dev.off()








#Figure 4G - FastMNN with all acute mice (A1-3)#######

cell.orders<-c("Proliferating CD4",
               "Th1 Proliferating CD4",
               "Cytolytic Proliferating CD4",
               "Th1/Cytolytic TCM",
               "Th1/Th17 TEM",
               "Th2",
               "Th9",
               "Tfh/Th22",
               "Tregs",
               "IL10RA+ Tregs")

cell.orders.clusters<-c("1","3","7","8","0","9","6","2","5","4")

seurat.list.act.MNN3$seurat_clusters_types<-case_when(
  seurat.list.act.MNN3$seurat_clusters == 1 ~ "Proliferating CD4",
  seurat.list.act.MNN3$seurat_clusters == 3 ~ "Th1 Proliferating CD4",
  seurat.list.act.MNN3$seurat_clusters == 7 ~ "Cytolytic Proliferating CD4",
  seurat.list.act.MNN3$seurat_clusters == 8 ~ "Th1/Cytolytic TCM",
  seurat.list.act.MNN3$seurat_clusters == 0 ~ "Th1/Th17 TEM",
  seurat.list.act.MNN3$seurat_clusters == 9 ~ "Th2",
  seurat.list.act.MNN3$seurat_clusters == 6 ~ "Th9",
  seurat.list.act.MNN3$seurat_clusters == 2 ~ "Tfh/Th22",
  seurat.list.act.MNN3$seurat_clusters == 5 ~ "Tregs",
  seurat.list.act.MNN3$seurat_clusters == 4 ~ "IL10RA+ Tregs",
)

#Let's extract the numbers of cells HIV-/HIV+ using cluster designations

total_vector<-seurat.list.act.MNN3@meta.data%>%
  group_by(seurat_clusters_types)%>%
  dplyr::count(status=="HIV-")

HIVpos.table<-total_vector[total_vector$`status == "HIV-"`=="FALSE",]
HIVneg.table<-total_vector[total_vector$`status == "HIV-"`=="TRUE",]

total_vector1<-seurat.list.act.MNN3@meta.data%>%
  group_by(seurat_clusters_types)%>%
  dplyr::count(seurat_clusters_types)

TableC <- data.frame(
  HIVpos=c(total_vector1$n-HIVneg.table$n),#Sucesses
  HIVneg=HIVneg.table$n#Failures
)

#Unfortunately I'm using a version of the seurat object which doesn't have the final cluster desigantions so I'll add them manually
TableC<-`rownames<-`(TableC,HIVneg.table$seurat_clusters_types)

##Lets do the Fisher's Exact Test

#Fishers Exact Test
TableC_fisher<-fisher_test(TableC,simulate.p.value = TRUE)

#Fishers Exact Test - Post hoc pairwise
pairwisetest_Fig3H<-pairwise_fisher_test(x=TableC,p.adjust.method="bonferroni")

#Comparing each cluster to the overall average

TableC_results<-tibble()
for(i in 1:nrow(TableC)) {
  TableC_Total<-data_frame(HIVpos=sum(TableC$HIVpos),HIVneg=sum(TableC$HIVneg))
  if (TableC$HIVpos[i]>1){
    print("Doing Fishers exact test")
    TableC1<-TableC[i,]
    TableC1<-rbind(TableC1,TableC_Total)
    res<-pairwise_fisher_test(TableC1,p.adjust.method = "bonferroni")
    print(res)
    TableC_results<-rbind(TableC_results,res)
  }
  else{
    print("No HIV+ -> skip")
  }
}

####Create two barplot versions
a<-ggplot(seurat.list.act.MNN3@meta.data,aes(x=seurat_clusters, fill=status))+
  geom_bar(position="fill")+
  geom_hline(yintercept = TableC_Total$HIVpos/(TableC_Total$HIVneg+TableC_Total$HIVpos), color="darkgrey",linetype=2)+
  scale_y_continuous(labels = scales::percent,expand=c(0,0))+
  scale_x_discrete(limits=cell.orders.clusters,labels=cell.orders)+
  theme_bw()+
  labs(fill="HIV status",y="% of cells",x="")+
  theme(axis.text.x = element_text(hjust=1,angle=45),legend.key.size = unit(0.1, 'cm'))+
  scale_fill_manual(values = HIV.cols,breaks = c("HIV+ high","HIV+ low","HIV-"))+
  coord_cartesian(ylim = c(0,0.1))

b<-ggplot(seurat.list.act.MNN3@meta.data,aes(x=seurat_clusters, fill=status))+
  geom_bar(stat ="count",position="stack")+
  scale_x_discrete(limits=cell.orders.clusters,labels=cell.orders)+
  theme_bw()+
  labs(fill="HIV status",y="# of cells",x="")+
  theme(axis.text.x = element_text(hjust=1,angle=45),legend.key.size = unit(0.1, 'cm'))+
  scale_fill_manual(values = HIV.cols,breaks = c("HIV+ high","HIV+ low","HIV-"))

####Create heatmaps to visualize pairwise results

c<-ggplot(pairwisetest_Fig3H,
          aes(group1, group2))+ # x and y axes => Var1 and Var2
  geom_tile(aes(fill = p.adj))+ # background colours are mapped according to the value column
  geom_text(aes(label = p.adj),size=2)+ # write the values
  scale_fill_gradient2(low = "darkred",
                       mid = "white",
                       high = "white",
                       midpoint=0.05,
                       limits=c(min(pairwisetest_Fig3H$p.adj),max(pairwisetest_Fig3H$p.adj)),
                       trans="log") + # determine the colour
  scale_x_discrete(limits=unique(pairwisetest_Fig3H$group1),name="")+
  scale_y_discrete(limits=unique(pairwisetest_Fig3H$group2),name="")+
  NoLegend()+
  theme(panel.grid.major.x=element_blank(), 
        panel.grid.minor.x=element_blank(), 
        panel.grid.major.y=element_blank(), 
        panel.grid.minor.y=element_blank(),
        panel.background=element_rect(fill="white"),
        axis.text.x = element_text(angle=45, hjust = 1,vjust=1,size = 8),
        axis.text.y = element_text(size = 8))

d<-ggplot(pairwisetest_Fig3H,
          aes(group1, group2))+ # x and y axes => Var1 and Var2
  geom_tile(aes(fill = p.adj.signif))+ # background colours are mapped according to the value column
  geom_text(aes(label = p.adj.signif))+ # write the values
  scale_fill_manual(values=c("****" = "darkred",
                             "***" = "red",
                             "**"="orange",
                             "*"="lightyellow",
                             "ns" = "white"))+ # determine the colour
  scale_x_discrete(limits=unique(pairwisetest_Fig3H$group1),name="")+
  scale_y_discrete(limits=unique(pairwisetest_Fig3H$group2),name="")+
  NoLegend()+
  theme(panel.grid.major.x=element_blank(), 
        panel.grid.minor.x=element_blank(), 
        panel.grid.major.y=element_blank(), 
        panel.grid.minor.y=element_blank(),
        panel.background=element_rect(fill="white"),
        axis.text.x = element_text(angle=45, hjust = 1,vjust=1,size = 8),
        axis.text.y = element_text(size = 8))

pdf("Statistics_Fig4G.pdf",height = 14)
ggarrange(a,b,ncol = 1,nrow = 2)
ggarrange(c,d,ncol = 1,nrow = 2)
dev.off()







#Figure 4I - Dataset with treated mice (T1/2)#######

seurat_treated2 <- FindClusters(seurat_treated,resolution = 0.5) 
seurat_treated2$seurat_clusters2<-case_when(seurat_treated2$seurat_clusters==0 ~ "Tfh-like CD4",
                                            seurat_treated2$seurat_clusters==1 ~ "Th1/Th17 Memory",
                                            seurat_treated2$seurat_clusters==2 ~ "Proliferating CD4",
                                            seurat_treated2$seurat_clusters==3 ~ "Cytolytic CD4 Memory",
                                            seurat_treated2$seurat_clusters==4 ~ "Naive CD4",
                                            seurat_treated2$seurat_clusters==5 ~ "Tfh/Th22",
                                            seurat_treated2$seurat_clusters==6 ~ "Tregs",
                                            seurat_treated2$seurat_clusters==7 ~ "MT Cells")


cell.orders2<-c("Proliferating CD4",
                "Naive CD4",
                "Cytolytic CD4 Memory",
                "Th1/Th17 Memory",
                "Tfh-like CD4",
                "Tfh/Th22",
                "Tregs",
                "MT Cells")

#Let's extract the numbers of cells HIV-/HIV+ using cluster designations

total_vector<-seurat_treated2@meta.data%>%
  group_by(seurat_clusters2)%>%
  dplyr::count(Fluorescence=="dsRed")

GFP.table<-total_vector[total_vector$`Fluorescence == "dsRed"`=="FALSE",]
dsRed.table<-total_vector[total_vector$`Fluorescence == "dsRed"`=="TRUE",]

total_vector1<-seurat_treated2@meta.data%>%
  group_by(seurat_clusters2)%>%
  dplyr::count(seurat_clusters2)

TableE <- data.frame(
  GFP=c(total_vector1$n-dsRed.table$n),#Sucesses
  dsRed=dsRed.table$n#Failures
)

#Unfortunately I'm using a version of the seurat object which doesn't have the final cluster desigantions so I'll add them manually
TableE<-`rownames<-`(TableE,dsRed.table$seurat_clusters2)

##Lets do the Fisher's Exact Test

#Fishers Exact Test
TableE_fisher<-fisher_test(TableE,simulate.p.value = TRUE)

#Fishers Exact Test - Post hoc pairwise
pairwisetest_Fig4I<-pairwise_fisher_test(x=TableE,p.adjust.method="bonferroni")

#Comparing each cluster to the overall average

TableE_results<-tibble()
for(i in 1:nrow(TableE)) {
  TableE_Total<-data_frame(GFP=sum(TableE$GFP),dsRed=sum(TableE$dsRed))
  if (TableE$GFP[i]>1){
    print("Doing Fishers exact test")
    TableE1<-TableE[i,]
    TableE1<-rbind(TableE1,TableE_Total)
    res<-pairwise_fisher_test(TableE1,p.adjust.method = "bonferroni")
    print(res)
    TableE_results<-rbind(TableE_results,res)
  }
  else{
    print("No GFP -> skip")
  }
}

####Create two barplot versions
a<-ggplot(seurat_treated2@meta.data,aes(x=seurat_clusters2, fill=Fluorescence))+
  geom_bar(position="fill")+
  geom_hline(yintercept = TableE_Total$GFP/(TableE_Total$dsRed+TableE_Total$GFP), color="darkgrey",linetype=2)+
  scale_y_continuous(labels = scales::percent,expand=c(0,0))+
  scale_x_discrete(limits=cell.orders2)+
  theme_bw()+
  labs(fill="Fluorescence",y="% of cells",x="")+
  theme(axis.text.x = element_text(hjust=1,angle=45),legend.key.size = unit(0.1, 'cm'))+
  scale_fill_manual(values = fluor.cols)+
  coord_cartesian(ylim = c(0,0.05))


b<-ggplot(seurat_treated2@meta.data,aes(x=seurat_clusters2, fill=Fluorescence))+
  geom_bar(stat = "count",position="stack")+
  scale_x_discrete(limits=cell.orders2)+
  theme_bw()+
  labs(fill="Fluorescence",y="# of cells",x="")+
  theme(axis.text.x = element_text(hjust=1,angle=45),legend.key.size = unit(0.1, 'cm'))+
  scale_fill_manual(values = fluor.cols)

####Create heatmaps to visualize pairwise results

c<-ggplot(pairwisetest_Fig4I,
          aes(group1, group2))+ # x and y axes => Var1 and Var2
  geom_tile(aes(fill = p.adj))+ # background colours are mapped according to the value column
  geom_text(aes(label = p.adj),size=2)+ # write the values
  scale_fill_gradient2(low = "darkred",
                       mid = "white",
                       high = "white",
                       midpoint=0.05,
                       limits=c(min(pairwisetest_Fig4I$p.adj),max(pairwisetest_Fig4I$p.adj)),
                       trans="log") + # determine the colour
  scale_x_discrete(limits=unique(pairwisetest_Fig4I$group1),name="")+
  scale_y_discrete(limits=unique(pairwisetest_Fig4I$group2),name="")+
  NoLegend()+
  theme(panel.grid.major.x=element_blank(), 
        panel.grid.minor.x=element_blank(), 
        panel.grid.major.y=element_blank(), 
        panel.grid.minor.y=element_blank(),
        panel.background=element_rect(fill="white"),
        axis.text.x = element_text(angle=45, hjust = 1,vjust=1,size = 8),
        axis.text.y = element_text(size = 8))

d<-ggplot(pairwisetest_Fig4I,
          aes(group1, group2))+ # x and y axes => Var1 and Var2
  geom_tile(aes(fill = p.adj.signif))+ # background colours are mapped according to the value column
  geom_text(aes(label = p.adj.signif))+ # write the values
  scale_fill_manual(values=c("****" = "darkred",
                             "***" = "red",
                             "**"="orange",
                             "*"="lightyellow",
                             "ns" = "white"))+ # determine the colour
  scale_x_discrete(limits=unique(pairwisetest_Fig4I$group1),name="")+
  scale_y_discrete(limits=unique(pairwisetest_Fig4I$group2),name="")+
  NoLegend()+
  theme(panel.grid.major.x=element_blank(), 
        panel.grid.minor.x=element_blank(), 
        panel.grid.major.y=element_blank(), 
        panel.grid.minor.y=element_blank(),
        panel.background=element_rect(fill="white"),
        axis.text.x = element_text(angle=45, hjust = 1,vjust=1,size = 8),
        axis.text.y = element_text(size = 8))

pdf("Statistics_Fig6E.pdf",height = 14)
ggarrange(a,b,ncol = 1,nrow = 2)
ggarrange(c,d,ncol = 1,nrow = 2)
dev.off()




