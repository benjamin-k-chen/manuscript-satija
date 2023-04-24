##Barplots using seurat.list.post.3rd.UMAP

seurat.list.post.3rdUMAP<-readRDS("~/Desktop/seurat.list.merged.by.mouse.copy.PostUMAP2023-04-20.rds")
seurat.list.post.3rdUMAP[[8]]$predicted.celltype.l2 <- str_replace_all(seurat.list.post.3rdUMAP[[8]]$predicted.celltype.l2,"CD4 Others Proliferating","CD4 Others")
seurat.list.post.3rdUMAP[[1]]$predicted.celltype.l2 <- str_replace_all(seurat.list.post.3rdUMAP[[1]]$predicted.celltype.l2,"CD4 Others Proliferating","CD4 Others")

#####
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
library(ggpubr)
#######
pdf("Barplots-20-04-23.pdf")
#BarPlots
for (i in 1:length(seurat.list.post.3rdUMAP)){
  
  #Barplots - predicted.celltype.l1 and l2
  
  #Creating a celltype.l2 table
  cell_table <- seurat.list.post.3rdUMAP[[i]]@meta.data%>%
    group_by(seurat_clusters)%>%
    dplyr::count(predicted.celltype.l2)%>%
    mutate(percent=(n/sum(n)))
  
  cell_table2<-cell_table[order(match(cell_table$predicted.celltype.l2, cell.order.rev)),]
  
  cell_table2$predicted.celltype.l2<-factor(cell_table2$predicted.celltype.l2,levels = unique(cell_table2$predicted.celltype.l2))
  
  a1<-ggplot(cell_table2,aes(x=seurat_clusters,y=cell_table2$percent,fill=factor(cell_table2$predicted.celltype.l2,levels = cell.order.rev)))+
    geom_bar(stat = "identity",position = "stack")+
    scale_y_continuous(expand=c(0,0),labels = scales::percent)+
    theme_bw()+labs(fill="predicted.celltype.l2",y="% of cells",x="Cluster",caption = paste0(dplyr::first(seurat.list.post.3rdUMAP[[i]]$Count),"-final barplots"))+
    theme(axis.text.x = element_text(hjust=1,angle=45),legend.key.size = unit(0.1, 'cm'))+
    scale_fill_manual(values = cell.cols,breaks = cell.order.rev)
  
  #cerating the .l1 barplots (-> celltypes are ordered alphabetically)
  a2<-ggplot(seurat.list.post.3rdUMAP[[i]]@meta.data,aes(x=seurat_clusters, fill=predicted.celltype.l1))+
    geom_bar(position="fill")+
    scale_y_continuous(labels = scales::percent,expand=c(0,0))+
    theme_bw()+
    labs(fill="predicted.celltype.l1",y="% of cells",x="Cluster",caption = paste0(dplyr::first(seurat.list.post.3rdUMAP[[i]]$Count),"-final barplots"))+
    theme(axis.text.x = element_text(hjust=1,angle=45),legend.key.size = unit(0.1, 'cm'))

  show(ggarrange(a1,a2,ncol=1,nrow = 2))
  
  #Barplots -HIV status
  
  h1<-ggplot(seurat.list.post.3rdUMAP[[i]]@meta.data,aes(x=seurat_clusters, fill=status)) + geom_bar(position="fill")+
    scale_y_continuous(labels = scales::percent,expand=c(0,0))+
    theme_bw()+
    labs(fill="HIV status",y="% of cells",x="Cluster",caption = paste0(dplyr::first(seurat.list.post.3rdUMAP[[i]]$Count),"-final barplots"))+
    theme(axis.text.x = element_text(hjust=1,angle=45),legend.key.size = unit(0.1, 'cm'))+
    scale_fill_manual(values = HIV.cols,breaks = c("HIV+ high","HIV+ low","HIV-"))+
    coord_cartesian(ylim = c(0,0.25))
  
  h2<-ggplot(seurat.list.post.3rdUMAP[[i]]@meta.data,aes(x=predicted.celltype.l2, fill=status)) + geom_bar(position="fill")+
    scale_y_continuous(labels = scales::percent,expand=c(0,0))+
    theme_bw()+
    labs(fill="HIV status",y="% of cells",x="predicted.celltypes.l2",caption = paste0(dplyr::first(seurat.list.post.3rdUMAP[[i]]$Count),"-final barplots"))+
    theme(axis.text.x = element_text(hjust=1,angle=45),legend.key.size = unit(0.1, 'cm'))+
    scale_fill_manual(values = HIV.cols,breaks = c("HIV+ high","HIV+ low","HIV-"))+
    coord_cartesian(ylim = c(0,0.25))
  show(ggarrange(h1,h2,ncol=1,nrow = 2))

  #Barplots - LentiRG expression
  
  t1<-ggplot(seurat.list.post.3rdUMAP[[i]]@meta.data,aes(x=seurat_clusters, fill=LentiRG.expression))+geom_bar(position="fill")+
    scale_y_continuous(labels = scales::percent,expand=c(0,0))+theme_bw()+
    labs(fill="LentiRG status",y="% of cells",x="Cluster",caption = paste0(dplyr::first(seurat.list.post.3rdUMAP[[i]]$Count),"-final barplots"))+
    theme(axis.text.x = element_text(hjust=1,angle=45),legend.key.size = unit(0.1, 'cm'))+
    scale_fill_manual(values = dsRed.cols)
  
  t2<-ggplot(seurat.list.post.3rdUMAP[[i]]@meta.data,aes(x=predicted.celltype.l2, fill=LentiRG.expression))+geom_bar(position="fill")+
    scale_y_continuous(labels = scales::percent,expand=c(0,0))+theme_bw()+
    labs(fill="LentiRG status",y="% of cells",x="predicted.celltypes.l2",caption = paste0(dplyr::first(seurat.list.post.3rdUMAP[[i]]$Count),"-final barplots"))+
    theme(axis.text.x = element_text(hjust=1,angle=45),legend.key.size = unit(0.1, 'cm'))+
    scale_fill_manual(values = dsRed.cols)
  
  show(ggarrange(t1,t2,ncol=1,nrow = 2))
  
  #Barplots - Fluorescence
  o1<-ggplot(seurat.list.post.3rdUMAP[[i]]@meta.data,aes(x=seurat_clusters, fill=Fluorescence))+geom_bar(position="fill")+
    scale_y_continuous(labels = scales::percent,expand=c(0,0))+theme_bw()+
    labs(fill="Fluorescence",y="% of cells",x="Cluster",caption = paste0(dplyr::first(seurat.list.post.3rdUMAP[[i]]$Count),"-final barplots"))+
    theme(axis.text.x = element_text(hjust=1,angle=45),legend.key.size = unit(0.1, 'cm'))+
    scale_fill_manual(values = fluor.cols)
  
  o2<-ggplot(seurat.list.post.3rdUMAP[[i]]@meta.data,aes(x=predicted.celltype.l2, fill=Fluorescence)) + geom_bar(position="fill")+
    scale_y_continuous(labels = scales::percent,expand=c(0,0))+theme_bw()+
    labs(fill="Fluorescence",y="% of cells",x="predicted.celltype.l2",caption = paste0(dplyr::first(seurat.list.post.3rdUMAP[[i]]$Count),"-final barplots"))+
    theme(axis.text.x = element_text(hjust=1,angle=45),legend.key.size = unit(0.1, 'cm'))+
    scale_fill_manual(values = fluor.cols)
  
  show(ggarrange(o1,o2,ncol=1,nrow = 2))
}
  
dev.off()




cell_table <- seurat.list.post.3rdUMAP[[3]]@meta.data%>%
  group_by(seurat_clusters)%>%
  dplyr::count(predicted.celltype.l2)%>%
  mutate(percent=(n/sum(n))*100)

cell_table2<-cell_table[order(match(cell_table$predicted.celltype.l2, cell.order.rev)),]

cell_table2$predicted.celltype.l2<-factor(cell_table2$predicted.celltype.l2,levels = unique(cell_table2$predicted.celltype.l2))



ggplot(seurat.list.post.3rdUMAP[[3]]@meta.data,aes(x=seurat_clusters, fill=predicted.celltype.l2))+
  geom_bar(position="fill")+
  scale_y_continuous(labels = scales::percent)+
  theme_bw()+
  labs(title = paste0(unique(seurat.list.post.3rdUMAP[[3]]$Count)))+
  theme(axis.text.x = element_text(hjust=1,angle=45),legend.key.size = unit(0.1, 'cm'))+
  scale_fill_manual(values = cell.cols,breaks = cell.order.rev)


ggplot(cell_table2,aes(x=seurat_clusters,y=cell_table2$percent,fill=factor(cell_table2$predicted.celltype.l2,levels = cell.order.rev)))+
  geom_bar(stat = "identity",position = "stack")+
  theme_bw()+labs(fill="predicted.celltype.l2")+
  theme(axis.text.x = element_text(hjust=1,angle=45),legend.key.size = unit(0.1, 'cm'))+
  scale_fill_manual(values = cell.cols,breaks = cell.order.rev,)






test.seurat<-seurat.list.post.3rdUMAP[[2]]

columnstoremove<-str_subset(colnames(seurat.list.post.3rdUMAP[[2]]@meta.data),pattern = "SCT_snn_res")
for(i in columnstoremove) {
  test.seurat[[i]] <- NULL
}

seurat<-lapply(X=seurat.list.post.3rdUMAP,FUN = SCTransform,vars.to.regress = "mitoHuCH38Percent")



DimPlot(seurat.list.post.3rdUMAP[[3]],dims = c(1,1))






pdf("20230407-Addiotional-3rdUMAP.pdf",width = 14,height = 10)
a1<-DimPlot(seurat.list.post.3rdUMAP[[3]],group.by = "Doublet.Singlet",pt.size = 0.5,order = c("Doublet","Singlet"),cols = c("Doublet"="darkblue","Singlet"="lightgrey"))

a2<-DimPlot(seurat.list.post.3rdUMAP[[7]],group.by = "Doublet.Singlet",pt.size = 0.5,order = c("Doublet","Singlet"),cols = c("Doublet"="darkblue","Singlet"="lightgrey"))

a3<-DimPlot(seurat.list.post.3rdUMAP[[8]],group.by = "Doublet.Singlet",pt.size = 0.5,order = c("Doublet","Singlet"),cols = c("Doublet"="darkblue","Singlet"="lightgrey"))

t1<-DimPlot(seurat.list.post.3rdUMAP[[4]],group.by = "Doublet.Singlet",pt.size = 0.5,order = c("Doublet","Singlet"),cols = c("Doublet"="darkblue","Singlet"="lightgrey"))

ggarrange(a1,a2,a3,t1,ncol =2, nrow = 2)
dev.off()









