#Setting the seed to allow reproducibility in terms of UMAPs
set.seed(10403)
#Defining the paths for the libraries use the packages from 
#Don't ask me why - the ones on native Minerva are sometimes outdated and using your personal one seem to be better (i.e. ggplot2)

#Answer these questions####
#1.What name do you want the output pdf to have?
pdf(file="TestSmallDE20230322.pdf",compress = TRUE)

#2.What is the filepath of the .RData workspace?
filePath <- c("/sc/arion/scratch/schmig04/Seurat_BC_test/Post3rdUMAP.RData")

#3.Where do you want the markers to be in?
OutPath <- c("/sc/arion/scratch/schmig04/Seurat_BC_test/markers20230322/")


#4.VolcanoPlot
Volc.cols<-c("grey","grey","grey","red")
Volc.labels = c("NS", "NS", "NS","fdr<0.05")

#Load Packages######
.libPaths(c("/hpc/packages/minerva-centos7/rpackages/4.2.0/site-library","/hpc/packages/minerva-centos7/rpackages/bioconductor/3.15","/hpc/users/schmig04/.Rlib","/hpc/packages/minerva-centos7/R/4.2.0/lib64/R/library"))

#
library(ggplot2,lib="/hpc/users/schmig04/.Rlib")
#
library(cli)
library(dplyr)
library(Seurat,lib="/hpc/packages/minerva-centos7/rpackages/4.2.0/site-library")
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


#Load the .RData workspace from the previous part#####
load(filePath)
Volc.labels = c("NS","NS","NS","fdr<0.05")

#Create Output directory#####
dir.create(OutPath)

#Differential gene expression analysis#######
for (i in 1:length(seurat.list.merged.by.mouse.copy.subset)){
  extracted.cells <- FetchData(object = seurat.list.merged.by.mouse.copy.subset[[i]], vars = c('orig.ident', "Count", 'nFeature_RNA', 'nCount_RNA',"Fluorescence","Treatment","SampleType","status"), slot="counts")
  extracted.cells$id <- rownames(extracted.cells)
  
  #GFP/dsRed/Unmarked cells
  
  if (c("dsRed") %in% unique(extracted.cells$Fluorescence) ){
    print(paste0(dplyr::first(seurat.list.merged.by.mouse.copy.subset[[i]]$Count)," - has dsRed cells"))
    dsRed.Cells <- extracted.cells[which(extracted.cells$Fluorescence == "dsRed"),]$id
    seurat.list.merged.by.mouse.copy.subset[[i]] <- SetIdent(object = seurat.list.merged.by.mouse.copy.subset[[i]], cells = dsRed.Cells  ,value = "dsRed")}
  
  if (c("Unmarked") %in% unique(extracted.cells$Fluorescence) ){
    print(paste0(dplyr::first(seurat.list.merged.by.mouse.copy.subset[[i]]$Count)," - has unmarked cells"))}
  
  if (c("GFP") %in% unique(extracted.cells$Fluorescence) ){
    print(paste0(dplyr::first(seurat.list.merged.by.mouse.copy.subset[[i]]$Count)," - has GFP cells"))
    GFP.Cells <- extracted.cells[which(extracted.cells$Fluorescence == "GFP"),]$id
    seurat.list.merged.by.mouse.copy.subset[[i]] <- SetIdent(object = seurat.list.merged.by.mouse.copy.subset[[i]], cells = GFP.Cells  ,value = "GFP")
    
    GFP.vs.dsRed.markers<-FindMarkers(seurat.list.merged.by.mouse.copy.subset[[i]], ident.1 = "dsRed", ident.2 = "GFP",test.use = "MAST", min.pct = 0.1, min.cells.group=1,logfc.threshold=0,  pseudocount.use= 0.001, random.seed = 12, latent.vars="nFeature_RNA")
    GFP.vs.dsRed.markers = GFP.vs.dsRed.markers[ !(rownames(GFP.vs.dsRed.markers) %in% featuresHIV), ]
    
    mousegenes<-rownames(GFP.vs.dsRed.markers)
    mousegenes<-str_subset(mousegenes,pattern = "mm10*")
    GFP.vs.dsRed.markers = GFP.vs.dsRed.markers[ !(rownames(GFP.vs.dsRed.markers) %in% mousegenes),]
    GFP.vs.dsRed.markers$fdr <- p.adjust(GFP.vs.dsRed.markers$p_val, method = "fdr", n = length(GFP.vs.dsRed.markers$p_val))
    
    show(EnhancedVolcano(GFP.vs.dsRed.markers, 
                         lab=rownames(GFP.vs.dsRed.markers),
                         x ="avg_log2FC", 
                         y ="p_val",subtitle = paste0("dsRed n=",nrow(extracted.cells[which(extracted.cells$Fluorescence %in% c("dsRed")),])," vs. GFP n=",nrow(extracted.cells[which(extracted.cells$Fluorescence %in% c("GFP")),])),
                         title = paste0("dsRed vs. GFP - ", dplyr::first(seurat.list.merged.by.mouse.copy.subset[[i]]$Count)),caption = "fdr cutoff = 0.05",
                         pCutoffCol = "fdr",pCutoff=0.05,FCcutoff = 0,col = Volc.cols,cutoffLineType = "blank",legendLabels = Volc.labels))
    
    write.csv(GFP.vs.dsRed.markers,file = paste0(OutPath,"dsRed-vs.-GFP-", dplyr::first(seurat.list.merged.by.mouse.copy.subset[[i]]$Count),"-MAST-markers.csv"))
    
    GFP.vs.dsRed.markers<-FindMarkers(seurat.list.merged.by.mouse.copy.subset[[i]], ident.1 = "dsRed", ident.2 = "GFP", min.pct = 0.1, min.cells.group=1,logfc.threshold=0,  pseudocount.use= 0.001, random.seed = 12, latent.vars="nFeature_RNA")
    GFP.vs.dsRed.markers = GFP.vs.dsRed.markers[ !(rownames(GFP.vs.dsRed.markers) %in% featuresHIV), ]
    
    mousegenes<-rownames(GFP.vs.dsRed.markers)
    mousegenes<-str_subset(mousegenes,pattern = "mm10*")
    GFP.vs.dsRed.markers = GFP.vs.dsRed.markers[ !(rownames(GFP.vs.dsRed.markers) %in% mousegenes),]
    GFP.vs.dsRed.markers$fdr <- p.adjust(GFP.vs.dsRed.markers$p_val, method = "fdr", n = length(GFP.vs.dsRed.markers$p_val))
    
    write.csv(GFP.vs.dsRed.markers,file = paste0(OutPath,"dsRed-vs.-GFP-", dplyr::first(seurat.list.merged.by.mouse.copy.subset[[i]]$Count),"-Wilcox-markers.csv"))
    
  }
  
  
  #HIV+/HIV-
  if (c("HIV-") %in% unique(extracted.cells$status) ){
    print(paste0(dplyr::first(seurat.list.merged.by.mouse.copy.subset[[i]]$Count)," - has HIV- cells"))
    HIVneg.Cells <- extracted.cells[which(extracted.cells$status == "HIV-"),]$id
    seurat.list.merged.by.mouse.copy.subset[[i]] <- SetIdent(object = seurat.list.merged.by.mouse.copy.subset[[i]], cells = HIVneg.Cells  ,value = "HIV-")
  }
  
  
  if (length(extracted.cells[which(extracted.cells$status %in% c("HIV+ very high","HIV+ high","HIV+ low")),]$status) > 1){
    print(paste0(dplyr::first(seurat.list.merged.by.mouse.copy.subset[[i]]$Count)," - has HIV+ cells"))
    HIVpos.Cells <- extracted.cells[which(extracted.cells$status %in% c("HIV+ very high","HIV+ high","HIV+ low")),]$id
    seurat.list.merged.by.mouse.copy.subset[[i]] <- SetIdent(object = seurat.list.merged.by.mouse.copy.subset[[i]], cells = HIVpos.Cells  ,value = "HIV+")
    
    
    HIVpos.vs.HIVneg.markers<-FindMarkers(seurat.list.merged.by.mouse.copy.subset[[i]], ident.1 = "HIV+", ident.2 = "HIV-",test.use = "MAST", min.pct = 0.1, min.cells.group=1,logfc.threshold=0,  pseudocount.use= 0.001, random.seed = 12, latent.vars="nFeature_RNA")
    HIVpos.vs.HIVneg.markers = HIVpos.vs.HIVneg.markers[ !(rownames(HIVpos.vs.HIVneg.markers) %in% featuresHIV), ]
    
    mousegenes<-rownames(HIVpos.vs.HIVneg.markers)
    mousegenes<-str_subset(mousegenes,pattern = "mm10*")
    HIVpos.vs.HIVneg.markers = HIVpos.vs.HIVneg.markers[ !(rownames(HIVpos.vs.HIVneg.markers) %in% mousegenes),]
    HIVpos.vs.HIVneg.markers$fdr <- p.adjust(HIVpos.vs.HIVneg.markers$p_val, method = "fdr", n = length(HIVpos.vs.HIVneg.markers$p_val))
    
    show(EnhancedVolcano(HIVpos.vs.HIVneg.markers, 
                         lab=rownames(HIVpos.vs.HIVneg.markers),
                         x ="avg_log2FC", 
                         y ="p_val",subtitle = paste0("HIV+ n=",nrow(extracted.cells[which(extracted.cells$status %in% c("HIV+ very high","HIV+ high","HIV+ low")),])," vs. HIV- n=",nrow(extracted.cells[which(extracted.cells$status %in% c("HIV-")),])),
                         title = paste0("HIV+ vs. HIV- - ", dplyr::first(seurat.list.merged.by.mouse.copy.subset[[i]]$Count)),caption = "fdr cutoff = 0.05",
                         pCutoffCol = "fdr",pCutoff=0.05,FCcutoff = 0,col = Volc.cols,cutoffLineType = "blank",legendLabels = Volc.labels))
    
    write.csv(HIVpos.vs.HIVneg.markers,file = paste0(OutPath,"HIVpos-vs.-HIVneg-", dplyr::first(seurat.list.merged.by.mouse.copy.subset[[i]]$Count),"-MAST-markers.csv"))
    
    HIVpos.vs.HIVneg.markers<-FindMarkers(seurat.list.merged.by.mouse.copy.subset[[i]], ident.1 = "HIV+", ident.2 = "HIV-", min.pct = 0.1, min.cells.group=1,logfc.threshold=0,  pseudocount.use= 0.001, random.seed = 12, latent.vars="nFeature_RNA")
    HIVpos.vs.HIVneg.markers = HIVpos.vs.HIVneg.markers[ !(rownames(HIVpos.vs.HIVneg.markers) %in% featuresHIV), ]
    
    mousegenes<-rownames(HIVpos.vs.HIVneg.markers)
    mousegenes<-str_subset(mousegenes,pattern = "mm10*")
    HIVpos.vs.HIVneg.markers = HIVpos.vs.HIVneg.markers[ !(rownames(HIVpos.vs.HIVneg.markers) %in% mousegenes),]
    HIVpos.vs.HIVneg.markers$fdr <- p.adjust(HIVpos.vs.HIVneg.markers$p_val, method = "fdr", n = length(HIVpos.vs.HIVneg.markers$p_val))
    
    write.csv(HIVpos.vs.HIVneg.markers,file = paste0(OutPath,"HIVpos-vs.-HIVneg-", dplyr::first(seurat.list.merged.by.mouse.copy.subset[[i]]$Count),"-Wilcox-markers.csv"))
    
    
  }
  
  
  if (c("HIV+ low") %in% unique(extracted.cells$status) ){print(paste0(dplyr::first(seurat.list.merged.by.mouse.copy.subset[[i]]$Count)," - has HIV+ low cells"))}
  
  if (c("HIV+ high") %in% unique(extracted.cells$status) ){
    print(paste0(dplyr::first(seurat.list.merged.by.mouse.copy.subset[[i]]$Count)," - has HIV+ high cells"))
    
    HIVpos.Cells <- extracted.cells[which(extracted.cells$status %in% c("HIV+ high")),]$id
    seurat.list.merged.by.mouse.copy.subset[[i]] <- SetIdent(object = seurat.list.merged.by.mouse.copy.subset[[i]], cells = HIVpos.Cells  ,value = "HIV+ high")
    
    
    HIVpos.vs.HIVneg.markers<-FindMarkers(seurat.list.merged.by.mouse.copy.subset[[i]], ident.1 = "HIV+ high", ident.2 = "HIV-",test.use = "MAST", min.pct = 0.1, min.cells.group=1,logfc.threshold=0,  pseudocount.use= 0.001, random.seed = 12, latent.vars="nFeature_RNA")
    HIVpos.vs.HIVneg.markers = HIVpos.vs.HIVneg.markers[ !(rownames(HIVpos.vs.HIVneg.markers) %in% featuresHIV), ]
    
    mousegenes<-rownames(HIVpos.vs.HIVneg.markers)
    mousegenes<-str_subset(mousegenes,pattern = "mm10*")
    HIVpos.vs.HIVneg.markers = HIVpos.vs.HIVneg.markers[ !(rownames(HIVpos.vs.HIVneg.markers) %in% mousegenes),]
    HIVpos.vs.HIVneg.markers$fdr <- p.adjust(HIVpos.vs.HIVneg.markers$p_val, method = "fdr", n = length(HIVpos.vs.HIVneg.markers$p_val))
    
    show(EnhancedVolcano(HIVpos.vs.HIVneg.markers, 
                         lab=rownames(HIVpos.vs.HIVneg.markers),
                         x ="avg_log2FC", 
                         y ="p_val",subtitle = paste0("HIV+ high n=",nrow(extracted.cells[which(extracted.cells$status %in% c("HIV+ high")),])," vs. HIV- n=",nrow(extracted.cells[which(extracted.cells$status %in% c("HIV-")),])),
                                                      title = paste0("HIV+ high vs. HIV- - ", dplyr::first(seurat.list.merged.by.mouse.copy.subset[[i]]$Count)),caption = "fdr cutoff = 0.05",
                                                      pCutoffCol = "fdr",pCutoff=0.05,FCcutoff = 0,col = Volc.cols,cutoffLineType = "blank",legendLabels = Volc.labels))
    write.csv(HIVpos.vs.HIVneg.markers,file = paste0(OutPath,"HIVhighpos-vs.-HIVneg-", dplyr::first(seurat.list.merged.by.mouse.copy.subset[[i]]$Count),"-MAST-markers.csv"))
    
    HIVpos.vs.HIVneg.markers<-FindMarkers(seurat.list.merged.by.mouse.copy.subset[[i]], ident.1 = "HIV+ high", ident.2 = "HIV-", min.pct = 0.1, min.cells.group=1,logfc.threshold=0,  pseudocount.use= 0.001, random.seed = 12, latent.vars="nFeature_RNA")
    HIVpos.vs.HIVneg.markers = HIVpos.vs.HIVneg.markers[ !(rownames(HIVpos.vs.HIVneg.markers) %in% featuresHIV), ]
    mousegenes<-rownames(HIVpos.vs.HIVneg.markers)
    mousegenes<-str_subset(mousegenes,pattern = "mm10*")
    HIVpos.vs.HIVneg.markers = HIVpos.vs.HIVneg.markers[ !(rownames(HIVpos.vs.HIVneg.markers) %in% mousegenes),]
    HIVpos.vs.HIVneg.markers$fdr <- p.adjust(HIVpos.vs.HIVneg.markers$p_val, method = "fdr", n = length(HIVpos.vs.HIVneg.markers$p_val))
    
    write.csv(HIVpos.vs.HIVneg.markers,file = paste0(OutPath,"HIVhighpos-vs.-HIVneg-", dplyr::first(seurat.list.merged.by.mouse.copy.subset[[i]]$Count),"-Wilcox-markers.csv"))
    }
  
  if (c("HIV+ very high") %in% unique(extracted.cells$status) ){
    print(paste0(dplyr::first(seurat.list.merged.by.mouse.copy.subset[[i]]$Count)," - has HIV+ very high cells"))
    
    HIVpos.Cells <- extracted.cells[which(extracted.cells$status %in% c("HIV+ very high")),]$id
    seurat.list.merged.by.mouse.copy.subset[[i]] <- SetIdent(object = seurat.list.merged.by.mouse.copy.subset[[i]], cells = HIVpos.Cells  ,value = "HIV+ very high")
    
    
    HIVpos.vs.HIVneg.markers<-FindMarkers(seurat.list.merged.by.mouse.copy.subset[[i]], ident.1 = "HIV+ very high", ident.2 = "HIV-",test.use = "MAST", min.pct = 0.1, min.cells.group=1,logfc.threshold=0,  pseudocount.use= 0.001, random.seed = 12, latent.vars="nFeature_RNA")
    HIVpos.vs.HIVneg.markers = HIVpos.vs.HIVneg.markers[ !(rownames(HIVpos.vs.HIVneg.markers) %in% featuresHIV), ]
    
    mousegenes<-rownames(HIVpos.vs.HIVneg.markers)
    mousegenes<-str_subset(mousegenes,pattern = "mm10*")
    HIVpos.vs.HIVneg.markers = HIVpos.vs.HIVneg.markers[ !(rownames(HIVpos.vs.HIVneg.markers) %in% mousegenes),]
    HIVpos.vs.HIVneg.markers$fdr <- p.adjust(HIVpos.vs.HIVneg.markers$p_val, method = "fdr", n = length(HIVpos.vs.HIVneg.markers$p_val))
    
    show(EnhancedVolcano(HIVpos.vs.HIVneg.markers, 
                         lab=rownames(HIVpos.vs.HIVneg.markers),
                         x ="avg_log2FC", 
                         y ="p_val",subtitle = paste0("HIV+ very high n=",nrow(extracted.cells[which(extracted.cells$status %in% c("HIV+ very high")),])," vs. HIV- n=",nrow(extracted.cells[which(extracted.cells$status %in% c("HIV-")),])),
                                                      title = paste0("HIV+ very high vs. HIV- - ", dplyr::first(seurat.list.merged.by.mouse.copy.subset[[i]]$Count)),caption = "fdr cutoff = 0.05",
                                                      pCutoffCol = "fdr",pCutoff=0.05,FCcutoff = 0,col = Volc.cols,cutoffLineType = "blank",legendLabels = Volc.labels))
         
    write.csv(HIVpos.vs.HIVneg.markers,file = paste0(OutPath,"HIVveryhighpos-vs.-HIVneg-", dplyr::first(seurat.list.merged.by.mouse.copy.subset[[i]]$Count),"-MAST-markers.csv"))
         
    HIVpos.vs.HIVneg.markers<-FindMarkers(seurat.list.merged.by.mouse.copy.subset[[i]], ident.1 = "HIV+ very high", ident.2 = "HIV-", min.pct = 0.1, min.cells.group=1,logfc.threshold=0,  pseudocount.use= 0.001, random.seed = 12, latent.vars="nFeature_RNA")
    HIVpos.vs.HIVneg.markers = HIVpos.vs.HIVneg.markers[ !(rownames(HIVpos.vs.HIVneg.markers) %in% featuresHIV), ]
         
    mousegenes<-rownames(HIVpos.vs.HIVneg.markers)
    mousegenes<-str_subset(mousegenes,pattern = "mm10*")
    HIVpos.vs.HIVneg.markers = HIVpos.vs.HIVneg.markers[ !(rownames(HIVpos.vs.HIVneg.markers) %in% mousegenes),]
    HIVpos.vs.HIVneg.markers$fdr <- p.adjust(HIVpos.vs.HIVneg.markers$p_val, method = "fdr", n = length(HIVpos.vs.HIVneg.markers$p_val))
    
    write.csv(HIVpos.vs.HIVneg.markers,file = paste0(OutPath,"HIVveryhighpos-vs.-HIVneg-", dplyr::first(seurat.list.merged.by.mouse.copy.subset[[i]]$Count),"-Wilcox-markers.csv"))
         
         
  }
  
}




#HIV umi counts plots#######
##HIV expression determination at the end
for(i in 1:length(seurat.list.merged.by.mouse.copy.subset)) {
  gene_HIV2 <- ifelse(str_detect(rownames(seurat.list.merged.by.mouse.copy.subset[[i]]@assays$RNA),"gag-pol|^pol|pol-vif-vpr-tat-rev|vpu-env|env|p2a-cre-ires-nef"),"HIV+", "HIV-")
  hivpos_inds2 <- gene_HIV2 == "HIV+"
  hivneg_inds2 <- gene_HIV2 == "HIV-"
  cell_hivstat2 <- tibble(n_hivpos_umi = Matrix::colSums(seurat.list.merged.by.mouse.copy.subset[[i]]@assays$RNA[hivpos_inds2,]),
                          n_hivneg_umi = Matrix::colSums(seurat.list.merged.by.mouse.copy.subset[[i]]@assays$RNA[hivneg_inds2,]),
                          tot_umi = Matrix::colSums(seurat.list.merged.by.mouse.copy.subset[[i]]@assays$RNA),
                          Status = seurat.list.merged.by.mouse.copy.subset[[i]]$status,
                          Count = seurat.list.merged.by.mouse.copy.subset[[i]]$Count,
                          Mean_hivpos_umi = base::mean(n_hivpos_umi))
  show(ggplot(cell_hivstat2, aes(n_hivpos_umi, n_hivneg_umi, color = Status,shape=Count)) +
         geom_point()+scale_x_log10()+scale_y_log10()+labs(title = paste0(unique(cell_hivstat2$Count),"- HIV+ vs. HIV- UMIs"))+scale_color_manual(values = HIV.cols)+theme_bw())
  
  assign(paste0(unique(seurat.list.merged.by.mouse.copy.subset[[i]]$Count),"-HIVumis_table"),cell_hivstat2)
  table<-cell_hivstat2%>%
    group_by(Count)%>%
    dplyr::count(Status)
  assign(paste0(unique(seurat.list.merged.by.mouse.copy.subset[[i]]$Count),"-HIVStatus_table"),table)
  
  
  table<-cell_hivstat2[cell_hivstat2$Status=="HIV+ high",]
  table<-tibble(Count=table$Count,
                Status=table$Status,
                Mean_hiv_high_umi = mean(table$n_hivpos_umi))
  assign(paste0(unique(seurat.list.merged.by.mouse.copy.subset[[i]]$Count),"-HIV+high_table"),table)
  
  table<-cell_hivstat2[cell_hivstat2$Status=="HIV+ low",]
  table<-tibble(Count=table$Count,
                Status=table$Status,
                Mean_hiv_low_umi = mean(table$n_hivpos_umi))
  assign(paste0(unique(seurat.list.merged.by.mouse.copy.subset[[i]]$Count),"-HIV+low_table"),table)
  
  cell_redstat2 <- tibble(Count = seurat.list.merged.by.mouse.copy.subset[[i]]$Count,
                          LentiRG_Status = seurat.list.merged.by.mouse.copy.subset[[i]]$LentiRG.expression,
                          Fluorescence = seurat.list.merged.by.mouse.copy.subset[[i]]$Fluorescence)
  table2<-cell_redstat2%>%
    group_by(Count)%>%
    dplyr::count(LentiRG_Status)
  
  table3<-cell_redstat2%>%
    group_by(Count)%>%
    dplyr::count(Fluorescence)
  assign(paste0(unique(seurat.list.merged.by.mouse.copy.subset[[i]]$Count),"-LentiRGStatus_table"),table2)
  assign(paste0(unique(seurat.list.merged.by.mouse.copy.subset[[i]]$Count),"-Fluorescence_table"),table3)
  
}

#Remove clutter from the workspace#AUTOMATE####

rm("609_610_count","611_612_count","a1","a2","cell_hivstat",
   "cell_species","celltypes.after.mm","celltypes.after.mm.split","dirs","dsRed.Cells","dsRed.Cells.cluster","extracted.cells","GFP.Cells",
   "GFP.Cells.cluster","h1","h2","HIneg.Cells.cluster","HIVneg_trans.Cells","HIVpos_trans.Cells","HIVpos.Cells.cluste","i","mouse.genes","mousegenes",
   "p1","p2","p3","p4","p5","p6","p7","p8","seurat.list","seurat.list.individual2","seurat.list.markers.subset","seurat.list.markers.top10.subset",
   "table","table1","table2","tableXYZ","tbl","tbl1","tbl2","tblxyz","test.features",
   "Unmarked.Cells","Unmarked.Cells.cluster","x")

#Save of the R workspace#####
save.image("TestSmallFull.RData")

#Finish pdf######
dev.off()
