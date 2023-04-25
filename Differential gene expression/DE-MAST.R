#R script used to perform DGE using MAST
#
#We used the .RData workspace from the main script which contains our datasets within a object.list of seurat objects
#This explains the automation of the analysis in a loop
#Additionally we also used the Wilcoxon test method, the default in Seurat, to compare with the MAST results

#Setting the seed to allow reproducibility in terms of UMAPs
set.seed(10403)

#Creating an Output.pdf file
pdf(file="DE-Output.pdf",compress = TRUE)

#Saving the filepath to the .RData file resulting from the main script
filePath <- c("PATH")

#Creating a directory to save the marker.csv files
OutPath <- c("PATH_II")


#Definig colors for the volcano plots
Volc.cols<-c("grey","grey","grey","red")
Volc.labels = c("NS", "NS", "NS","fdr<0.05")

#Load Packages##
library(ggplot2)
library(cli)
library(dplyr)
library(Seurat) #We used Seurat v.4.1.0 since newere version are not compatible with MAST and no bugfix has been implmented yet
library(patchwork)
library(tidyverse)
library(ggrepel)
library(gridExtra)
library(ggpmisc)
library(remotes)
library(SeuratDisk)
library(ggpubr)
library(EnhancedVolcano)
library(MAST)


#Load the .RData workspace from the previous part#
load(filePath)

#Create Output directory##
dir.create(OutPath)

#Differential gene expression analysis#######
for (i in 1:length(seurat.list)){
  extracted.cells <- FetchData(object = seurat.list[[i]], vars = c('orig.ident', "Count", 'nFeature_RNA', 'nCount_RNA',"Fluorescence","Treatment","SampleType","status"), slot="counts")
  extracted.cells$id <- rownames(extracted.cells)
  
  #GFP/dsRed/Unmarked cells
  
  if (c("dsRed") %in% unique(extracted.cells$Fluorescence) ){
    print(paste0(dplyr::first(seurat.list[[i]]$Count)," - has dsRed cells"))
    dsRed.Cells <- extracted.cells[which(extracted.cells$Fluorescence == "dsRed"),]$id
    seurat.list[[i]] <- SetIdent(object = seurat.list[[i]], cells = dsRed.Cells  ,value = "dsRed")}
  
  if (c("Unmarked") %in% unique(extracted.cells$Fluorescence) ){
    print(paste0(dplyr::first(seurat.list[[i]]$Count)," - has unmarked cells"))}
  
  if (c("GFP") %in% unique(extracted.cells$Fluorescence) ){
    print(paste0(dplyr::first(seurat.list[[i]]$Count)," - has GFP cells"))
    GFP.Cells <- extracted.cells[which(extracted.cells$Fluorescence == "GFP"),]$id
    seurat.list[[i]] <- SetIdent(object = seurat.list[[i]], cells = GFP.Cells  ,value = "GFP")
    
    GFP.vs.dsRed.markers<-FindMarkers(seurat.list[[i]], ident.1 = "dsRed", ident.2 = "GFP",test.use = "MAST", min.pct = 0.1, min.cells.group=1,logfc.threshold=0,  pseudocount.use= 0.001, random.seed = 12, latent.vars="nFeature_RNA")
    GFP.vs.dsRed.markers = GFP.vs.dsRed.markers[ !(rownames(GFP.vs.dsRed.markers) %in% featuresHIV), ]
    
    mousegenes<-rownames(GFP.vs.dsRed.markers)
    mousegenes<-str_subset(mousegenes,pattern = "mm10*")
    GFP.vs.dsRed.markers = GFP.vs.dsRed.markers[ !(rownames(GFP.vs.dsRed.markers) %in% mousegenes),]
    GFP.vs.dsRed.markers$fdr <- p.adjust(GFP.vs.dsRed.markers$p_val, method = "fdr", n = length(GFP.vs.dsRed.markers$p_val))
    
    show(EnhancedVolcano(GFP.vs.dsRed.markers, 
                         lab=rownames(GFP.vs.dsRed.markers),
                         x ="avg_log2FC", 
                         y ="p_val",subtitle = paste0("dsRed n=",nrow(extracted.cells[which(extracted.cells$Fluorescence %in% c("dsRed")),])," vs. GFP n=",nrow(extracted.cells[which(extracted.cells$Fluorescence %in% c("GFP")),])),
                         title = paste0("dsRed vs. GFP - ", dplyr::first(seurat.list[[i]]$Count)),caption = "fdr cutoff = 0.05",
                         pCutoffCol = "fdr",pCutoff=0.05,FCcutoff = 0,col = Volc.cols,cutoffLineType = "blank",legendLabels = Volc.labels))
    
    write.csv(GFP.vs.dsRed.markers,file = paste0(OutPath,"dsRed-vs.-GFP-", dplyr::first(seurat.list[[i]]$Count),"-MAST-markers.csv"))
    
    GFP.vs.dsRed.markers<-FindMarkers(seurat.list[[i]], ident.1 = "dsRed", ident.2 = "GFP", min.pct = 0.1, min.cells.group=1,logfc.threshold=0,  pseudocount.use= 0.001, random.seed = 12, latent.vars="nFeature_RNA")
    GFP.vs.dsRed.markers = GFP.vs.dsRed.markers[ !(rownames(GFP.vs.dsRed.markers) %in% featuresHIV), ]
    
    mousegenes<-rownames(GFP.vs.dsRed.markers)
    mousegenes<-str_subset(mousegenes,pattern = "mm10*")
    GFP.vs.dsRed.markers = GFP.vs.dsRed.markers[ !(rownames(GFP.vs.dsRed.markers) %in% mousegenes),]
    GFP.vs.dsRed.markers$fdr <- p.adjust(GFP.vs.dsRed.markers$p_val, method = "fdr", n = length(GFP.vs.dsRed.markers$p_val))
    
    write.csv(GFP.vs.dsRed.markers,file = paste0(OutPath,"dsRed-vs.-GFP-", dplyr::first(seurat.list[[i]]$Count),"-Wilcox-markers.csv"))
    
  }
  
  
  #HIV+/HIV-
  if (c("HIV-") %in% unique(extracted.cells$status) ){
    print(paste0(dplyr::first(seurat.list[[i]]$Count)," - has HIV- cells"))
    HIVneg.Cells <- extracted.cells[which(extracted.cells$status == "HIV-"),]$id
    seurat.list[[i]] <- SetIdent(object = seurat.list[[i]], cells = HIVneg.Cells  ,value = "HIV-")
  }
  
  
  if (length(extracted.cells[which(extracted.cells$status %in% c("HIV+ high","HIV+ low")),]$status) > 1){
    print(paste0(dplyr::first(seurat.list[[i]]$Count)," - has HIV+ cells"))
    HIVpos.Cells <- extracted.cells[which(extracted.cells$status %in% c("HIV+ high","HIV+ low")),]$id
    seurat.list[[i]] <- SetIdent(object = seurat.list[[i]], cells = HIVpos.Cells  ,value = "HIV+")
    
    
    HIVpos.vs.HIVneg.markers<-FindMarkers(seurat.list[[i]], ident.1 = "HIV+", ident.2 = "HIV-",test.use = "MAST", min.pct = 0.1, min.cells.group=1,logfc.threshold=0,  pseudocount.use= 0.001, random.seed = 12, latent.vars="nFeature_RNA")
    HIVpos.vs.HIVneg.markers = HIVpos.vs.HIVneg.markers[ !(rownames(HIVpos.vs.HIVneg.markers) %in% featuresHIV), ]
    
    mousegenes<-rownames(HIVpos.vs.HIVneg.markers)
    mousegenes<-str_subset(mousegenes,pattern = "mm10*")
    HIVpos.vs.HIVneg.markers = HIVpos.vs.HIVneg.markers[ !(rownames(HIVpos.vs.HIVneg.markers) %in% mousegenes),]
    HIVpos.vs.HIVneg.markers$fdr <- p.adjust(HIVpos.vs.HIVneg.markers$p_val, method = "fdr", n = length(HIVpos.vs.HIVneg.markers$p_val))
    
    show(EnhancedVolcano(HIVpos.vs.HIVneg.markers, 
                         lab=rownames(HIVpos.vs.HIVneg.markers),
                         x ="avg_log2FC", 
                         y ="p_val",subtitle = paste0("HIV+ n=",nrow(extracted.cells[which(extracted.cells$status %in% c("HIV+ high","HIV+ low")),])," vs. HIV- n=",nrow(extracted.cells[which(extracted.cells$status %in% c("HIV-")),])),
                         title = paste0("HIV+ vs. HIV- - ", dplyr::first(seurat.list[[i]]$Count)),caption = "fdr cutoff = 0.05",
                         pCutoffCol = "fdr",pCutoff=0.05,FCcutoff = 0,col = Volc.cols,cutoffLineType = "blank",legendLabels = Volc.labels))
    
    write.csv(HIVpos.vs.HIVneg.markers,file = paste0(OutPath,"HIVpos-vs.-HIVneg-", dplyr::first(seurat.list[[i]]$Count),"-MAST-markers.csv"))
    
    HIVpos.vs.HIVneg.markers<-FindMarkers(seurat.list[[i]], ident.1 = "HIV+", ident.2 = "HIV-", min.pct = 0.1, min.cells.group=1,logfc.threshold=0,  pseudocount.use= 0.001, random.seed = 12, latent.vars="nFeature_RNA")
    HIVpos.vs.HIVneg.markers = HIVpos.vs.HIVneg.markers[ !(rownames(HIVpos.vs.HIVneg.markers) %in% featuresHIV), ]
    
    mousegenes<-rownames(HIVpos.vs.HIVneg.markers)
    mousegenes<-str_subset(mousegenes,pattern = "mm10*")
    HIVpos.vs.HIVneg.markers = HIVpos.vs.HIVneg.markers[ !(rownames(HIVpos.vs.HIVneg.markers) %in% mousegenes),]
    HIVpos.vs.HIVneg.markers$fdr <- p.adjust(HIVpos.vs.HIVneg.markers$p_val, method = "fdr", n = length(HIVpos.vs.HIVneg.markers$p_val))
    
    write.csv(HIVpos.vs.HIVneg.markers,file = paste0(OutPath,"HIVpos-vs.-HIVneg-", dplyr::first(seurat.list[[i]]$Count),"-Wilcox-markers.csv"))
    
    
  }
  
  
  if (c("HIV+ low") %in% unique(extracted.cells$status) ){print(paste0(dplyr::first(seurat.list[[i]]$Count)," - has HIV+ low cells"))}
  
  if (c("HIV+ high") %in% unique(extracted.cells$status) ){
    print(paste0(dplyr::first(seurat.list[[i]]$Count)," - has HIV+ high cells"))
    
    HIVpos.Cells <- extracted.cells[which(extracted.cells$status %in% c("HIV+ high")),]$id
    seurat.list[[i]] <- SetIdent(object = seurat.list[[i]], cells = HIVpos.Cells  ,value = "HIV+ high")
    
    
    HIVpos.vs.HIVneg.markers<-FindMarkers(seurat.list[[i]], ident.1 = "HIV+ high", ident.2 = "HIV-",test.use = "MAST", min.pct = 0.1, min.cells.group=1,logfc.threshold=0,  pseudocount.use= 0.001, random.seed = 12, latent.vars="nFeature_RNA")
    HIVpos.vs.HIVneg.markers = HIVpos.vs.HIVneg.markers[ !(rownames(HIVpos.vs.HIVneg.markers) %in% featuresHIV), ]
    
    mousegenes<-rownames(HIVpos.vs.HIVneg.markers)
    mousegenes<-str_subset(mousegenes,pattern = "mm10*")
    HIVpos.vs.HIVneg.markers = HIVpos.vs.HIVneg.markers[ !(rownames(HIVpos.vs.HIVneg.markers) %in% mousegenes),]
    HIVpos.vs.HIVneg.markers$fdr <- p.adjust(HIVpos.vs.HIVneg.markers$p_val, method = "fdr", n = length(HIVpos.vs.HIVneg.markers$p_val))
    
    show(EnhancedVolcano(HIVpos.vs.HIVneg.markers, 
                         lab=rownames(HIVpos.vs.HIVneg.markers),
                         x ="avg_log2FC", 
                         y ="p_val",subtitle = paste0("HIV+ high n=",nrow(extracted.cells[which(extracted.cells$status %in% c("HIV+ high")),])," vs. HIV- n=",nrow(extracted.cells[which(extracted.cells$status %in% c("HIV-")),])),
                                                      title = paste0("HIV+ high vs. HIV- - ", dplyr::first(seurat.list[[i]]$Count)),caption = "fdr cutoff = 0.05",
                                                      pCutoffCol = "fdr",pCutoff=0.05,FCcutoff = 0,col = Volc.cols,cutoffLineType = "blank",legendLabels = Volc.labels))
    write.csv(HIVpos.vs.HIVneg.markers,file = paste0(OutPath,"HIVhighpos-vs.-HIVneg-", dplyr::first(seurat.list[[i]]$Count),"-MAST-markers.csv"))
    
    HIVpos.vs.HIVneg.markers<-FindMarkers(seurat.list[[i]], ident.1 = "HIV+ high", ident.2 = "HIV-", min.pct = 0.1, min.cells.group=1,logfc.threshold=0,  pseudocount.use= 0.001, random.seed = 12, latent.vars="nFeature_RNA")
    HIVpos.vs.HIVneg.markers = HIVpos.vs.HIVneg.markers[ !(rownames(HIVpos.vs.HIVneg.markers) %in% featuresHIV), ]
    mousegenes<-rownames(HIVpos.vs.HIVneg.markers)
    mousegenes<-str_subset(mousegenes,pattern = "mm10*")
    HIVpos.vs.HIVneg.markers = HIVpos.vs.HIVneg.markers[ !(rownames(HIVpos.vs.HIVneg.markers) %in% mousegenes),]
    HIVpos.vs.HIVneg.markers$fdr <- p.adjust(HIVpos.vs.HIVneg.markers$p_val, method = "fdr", n = length(HIVpos.vs.HIVneg.markers$p_val))
    
    write.csv(HIVpos.vs.HIVneg.markers,file = paste0(OutPath,"HIVhighpos-vs.-HIVneg-", dplyr::first(seurat.list[[i]]$Count),"-Wilcox-markers.csv"))
    }
  
}


#HIV umi counts plots###
for(i in 1:length(seurat.list)) {
  gene_HIV2 <- ifelse(str_detect(rownames(seurat.list[[i]]@assays$RNA),"gag-pol|^pol|pol-vif-vpr-tat-rev|vpu-env|env|p2a-cre-ires-nef"),"HIV+", "HIV-")
  hivpos_inds2 <- gene_HIV2 == "HIV+"
  hivneg_inds2 <- gene_HIV2 == "HIV-"
  cell_hivstat2 <- tibble(n_hivpos_umi = Matrix::colSums(seurat.list[[i]]@assays$RNA[hivpos_inds2,]),
                          n_hivneg_umi = Matrix::colSums(seurat.list[[i]]@assays$RNA[hivneg_inds2,]),
                          tot_umi = Matrix::colSums(seurat.list[[i]]@assays$RNA),
                          Status = seurat.list[[i]]$status,
                          Count = seurat.list[[i]]$Count,
                          Mean_hivpos_umi = base::mean(n_hivpos_umi))
  show(ggplot(cell_hivstat2, aes(n_hivpos_umi, n_hivneg_umi, color = Status,shape=Count)) +
         geom_point()+scale_x_log10()+scale_y_log10()+labs(title = paste0(unique(cell_hivstat2$Count),"- HIV+ vs. HIV- UMIs"))+scale_color_manual(values = HIV.cols)+theme_bw())
  
  assign(paste0(unique(seurat.list[[i]]$Count),"-HIVumis_table"),cell_hivstat2)
  table<-cell_hivstat2%>%
    group_by(Count)%>%
    dplyr::count(Status)
  assign(paste0(unique(seurat.list[[i]]$Count),"-HIVStatus_table"),table)
  
  
  table<-cell_hivstat2[cell_hivstat2$Status=="HIV+ high",]
  table<-tibble(Count=table$Count,
                Status=table$Status,
                Mean_hiv_high_umi = mean(table$n_hivpos_umi))
  assign(paste0(unique(seurat.list[[i]]$Count),"-HIV+high_table"),table)
  
  table<-cell_hivstat2[cell_hivstat2$Status=="HIV+ low",]
  table<-tibble(Count=table$Count,
                Status=table$Status,
                Mean_hiv_low_umi = mean(table$n_hivpos_umi))
  assign(paste0(unique(seurat.list[[i]]$Count),"-HIV+low_table"),table)
  
  cell_redstat2 <- tibble(Count = seurat.list[[i]]$Count,
                          LentiRG_Status = seurat.list[[i]]$LentiRG.expression,
                          Fluorescence = seurat.list[[i]]$Fluorescence)
  table2<-cell_redstat2%>%
    group_by(Count)%>%
    dplyr::count(LentiRG_Status)
  
  table3<-cell_redstat2%>%
    group_by(Count)%>%
    dplyr::count(Fluorescence)
  assign(paste0(unique(seurat.list[[i]]$Count),"-LentiRGStatus_table"),table2)
  assign(paste0(unique(seurat.list[[i]]$Count),"-Fluorescence_table"),table3)
  
}


#Save of the R workspace#####
save.image("MAST-DE.RData")

#Finish pdf###
dev.off()


##END OF SCRIPT####
