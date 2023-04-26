#Seurat downstream analysis - BASIC WORKFLOW
#
#
#This is the basic workflow for the analysis of the paper Satija et al. (2023)
#
#
#
#Setting the seed to allow reproducibility in terms of UMAPs
set.seed(10403)

#Load packages we need####
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(DoubletFinder)
library(clustree)
library(scCustomize)

#Loading the dataset####

#Get filepath
filepath<-c('path/to/datasets')

dirs <- list.dirs(path = filepath, recursive = F, full.names = F)
dirs <- dirs[dirs %in% datasets]

#Loading all other datasets and creating corresponding SeuratObjects
#The name of the object is added as meta.data column 'Count'

for(x in dirs){
  name <-gsub('/filtered_feature_bc_matrix/','', x)
  name2 <-paste0(refnumber,name,sep="x")
  cts <- ReadMtx(mtx = paste0(filepath,x,'/matrix.mtx.gz'),
                 features = paste0(filepath,x,'/features.tsv.gz'),
                 cells = paste0(filepath,x,'/barcodes.tsv.gz'))
  
  # create seurat objects
  SeuratObject<-CreateSeuratObject(counts = cts)
  SeuratObject$Count <- paste0(name)
  assign(name2,SeuratObject)
}

rm(name,name2,cts,SeuratObject,x)

#Moving all SeuratObjects into a list of Objects
list <-ls(pattern="Ref5x")
seurat.list<-lapply(list,get)

#And remove the individual datasets from the enviroment
remove(list=ls(pattern = refnumber))

#Mitochondrial gene percentage#####
##Our custom reference was was made using the 10x reference that contains both human and mouse genes,
#therefore we can calculate both mitochondrial gene percentages

for(i in 1:length(seurat.list)) {
  seurat.list[[i]]$mitoHuCH38Percent <- PercentageFeatureSet(seurat.list[[i]], pattern='*-MT-') #Calculating human mitochondrial gene percentage
  seurat.list[[i]]$mitoMM10Percent <- PercentageFeatureSet(seurat.list[[i]], pattern='*-mt-') #Calculating mouse mitochondrial gene percentage
  
  p1<-VlnPlot(seurat.list[[i]], features=c("nCount_RNA","mitoMM10Percent","mitoHuCH38Percent","nFeature_RNA"),ncol=4,group.by = "Count")+labs(title = paste0(unique(seurat.list[[i]]$Count),"-before-QC"))
  p2<-FeatureScatter(seurat.list[[i]], feature1 = "nCount_RNA", feature2 = "mitoHuCH38Percent",group.by = "Count")+labs(title = paste0(unique(seurat.list[[i]]$Count),"-before-QC"))
  p3<-FeatureScatter(seurat.list[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "Count")+labs(title = paste0(unique(seurat.list[[i]]$Count),"-before-QC"))
  show(ggarrange(p1,p2,p3,ncol=2,nrow = 2))
}

#The individual threshold values can be found in the supplemetary data

#Metadata######
##Add recorded meta.data to our dataset
#We added Fluorescence-SampleType-Count-MouseID-Condition-Treatment
#The individual values can be found in the supplementary data

#Add Treatment as metadata

datasetsTRtbl<-data.frame(datasetsTR)
datasetsTRtbl$Count<-rownames(datasetsTRtbl)

for(i in 1:length(seurat.list)) {
  seurat.info <- tibble(Count = seurat.list[[i]]$Count)
  seurat.info <- inner_join(seurat.info,datasetsTRtbl,by="Count") 
  seurat.list[[i]] <- AddMetaData(seurat.list[[i]], metadata = seurat.info$datasetsTR, col.name = "Treatment")
}

#Fluorescence is added as metadata

datasetsFLtbl<-data.frame(datasetsFL)
datasetsFLtbl$Count<-rownames(datasetsFLtbl)

for(i in 1:length(seurat.list)) {
  seurat.info <- tibble(Count = seurat.list[[i]]$Count)
  seurat.info <- inner_join(seurat.info,datasetsFLtbl,by="Count") 
  seurat.list[[i]] <- AddMetaData(seurat.list[[i]], metadata = seurat.info$datasetsFL, col.name = "Fluorescence")
}


#MouseID is added as metadata

datasetsMouseIDtbl<-data.frame(datasetsMouseID)
datasetsMouseIDtbl$Count<-rownames(datasetsMouseIDtbl)

for(i in 1:length(seurat.list)) {
  seurat.info <- tibble(Count = seurat.list[[i]]$Count)
  seurat.info <- inner_join(seurat.info,datasetsMouseIDtbl,by="Count") 
  seurat.list[[i]] <- AddMetaData(seurat.list[[i]], metadata = seurat.info$datasetsFL, col.name = "MouseID")
}

#Add SampleType as metadata

datasetsSTtbl<-data.frame(datasetsST)
datasetsSTtbl$Count<-rownames(datasetsSTtbl)

for(i in 1:length(seurat.list)) {
  seurat.info <- tibble(Count = seurat.list[[i]]$Count)
  seurat.info <- inner_join(seurat.info,datasetsSTtbl,by="Count") 
  seurat.list[[i]] <- AddMetaData(seurat.list[[i]], metadata = seurat.info$datasetsST, col.name = "SampleType")
}

#Add Condition as metadata
datasetsCotbl<-data.frame(datasetsCo)
datasetsCotbl$Count<-rownames(datasetsCotbl)

for(i in 1:length(seurat.list)) {
  seurat.info <- tibble(Count = seurat.list[[i]]$Count)
  seurat.info <- inner_join(seurat.info,datasetsCotbl,by="Count") 
  seurat.list[[i]] <- AddMetaData(seurat.list[[i]], metadata = seurat.info$datasetsCo, col.name = "Condition")
}


#Species Identification#### 
#If you used a the 10X reference that contains both mouse and human genes

for(i in 1:length(seurat.list)){
  #We define every gene that has the prefix "mm10" as "mouse" and every other as "human"
  gene_species <- ifelse(str_detect(rownames(seurat.list[[i]]),"mm10*"), "mouse", "human")
  mouse_inds <- gene_species == "mouse"
  human_inds <- gene_species == "human"
  
  #Using these indicators we can now create a table calles "cell_species" using tibble
  #We calculate the amount of "human" or "mouse" umi based on the amount of genes that are "TRUE" within the total umi count
  #of each cell which we calculate using colSums() 
  #Within this table we also calculate the % of human/mouse umi within each cell (each row represents a cell)
  
  cell_species <- tibble(n_mouse_umi = Matrix::colSums(seurat.list[[i]][mouse_inds,]),
                         n_human_umi = Matrix::colSums(seurat.list[[i]][human_inds,]),
                         tot_umi = Matrix::colSums(seurat.list[[i]]),
                         prop_mouse = n_mouse_umi / tot_umi,
                         prop_human = n_human_umi / tot_umi,
                         Count = seurat.list[[i]]$Count)
  
  #Using our table we now determine the species of the cell in a new added column called "species"
  #With case_when() we decide that if prop_mouse bigger that >95% the value in species is to be "mouse"
  #if prop_human > 95% the value in species is to be "human", and if it's somewhere inbetween (e.g. 80% human, 20% mouse) the value is "mixed"
  
  cell_species <- cell_species %>% mutate(species = case_when(prop_mouse > 0.95 ~ "mouse",prop_human > 0.95 ~ "human",TRUE ~ "mixed"))
  table1<-cell_species %>% dplyr::count(species) %>% mutate(proportion = n / ncol(seurat.list[[i]]))
  show(ggplot(cell_species, aes(n_human_umi, n_mouse_umi, color = species,shape=Count)) +
         geom_point()+labs(title = paste0(unique(cell_species$Count),"- Human vs. Mouse UMIs"))+theme_bw())
  #We then add this information as metadata to our object
  seurat.list[[i]] <- AddMetaData(seurat.list[[i]], metadata = cell_species$species, col.name = "species")
  tbl1<-tableGrob(table1,rows=NULL,theme=tt)
  grid.arrange(ggplot(seurat.list[[i]]@meta.data,aes(x=Count, fill=species)) + geom_bar(position="fill")+
                 scale_y_continuous(labels = scales::percent)+theme_bw()+theme(axis.text.x = element_text(hjust=1,angle=45)), tbl1,ncol=2,as.table=TRUE,heights=c(2,1))
}




#HIV transcriptional positive cell identification####
#This step basically follows the same principle as the species identification step
#We initially defined a third category of 'HIV+' called 'HIV+ very high' this was relabeled as 'HIV+ high'

for(i in 1:length(seurat.list)) {
  gene_HIV <- ifelse(str_detect(rownames(seurat.list[[i]]),"mm10*|GRCh38*|dsRed|EGFP|WPRE"),"HIV-", "HIV+")
  hivpos_inds <- gene_HIV == "HIV+"
  hivneg_inds <- gene_HIV == "HIV-"
  cell_hivstat <- tibble(n_hivpos_umi = Matrix::colSums(seurat.list[[i]][hivpos_inds,]),
                         n_hivneg_umi = Matrix::colSums(seurat.list[[i]][hivneg_inds,]),
                         tot_umi = Matrix::colSums(seurat.list[[i]]),
                         prop_hivpos = n_hivpos_umi / tot_umi,
                         prop_hivneg = n_hivneg_umi / tot_umi,
                         SampleType = seurat.list[[i]]$SampleType,
                         Count = seurat.list[[i]]$Count)
  cell_hivstat <- cell_hivstat %>% mutate(status = case_when(n_hivpos_umi > 100 ~ "HIV+ very high",n_hivpos_umi > 10 & n_hivpos_umi <= 100 ~ "HIV+ high",n_hivpos_umi <= 1 ~ "HIV-",n_hivpos_umi > 1 & n_hivpos_umi <= 10 ~ "HIV+ low"))
  table2<-cell_hivstat %>% dplyr::count(status) %>% mutate(proportion = n / ncol(seurat.list[[i]]))
  show(ggplot(cell_hivstat, aes(n_hivpos_umi, n_hivneg_umi, color = status,shape=Count)) +
         geom_point()+scale_x_log10()+scale_y_log10()+labs(title = paste0(unique(cell_hivstat$Count),"- HIV+ vs. HIV- UMIs"))+scale_color_manual(values = HIV.cols)+theme_bw())
  seurat.list[[i]] <- AddMetaData(seurat.list[[i]], metadata = cell_hivstat$status, col.name = "status")
  tbl2<-tableGrob(table2,rows=NULL,theme=tt)
  grid.arrange(ggplot(seurat.list[[i]]@meta.data,aes(x=Count, fill=status)) + geom_bar(position="fill")+
                 scale_y_continuous(labels = scales::percent)+theme_bw()+theme(axis.text.x = element_text(hjust=1,angle=45)), tbl2,ncol=2,as.table=TRUE,heights=c(2,1))
}



#Lenti-RG transcriptionally positive cell identification####
#This follows the same principal as the two previous steps

for(i in 1:length(seurat.list)) {
  gene_dsRed <- ifelse(str_detect(rownames(seurat.list[[i]]),c("dsRed|WPRE|EGFP")),"LentiRG", "other")
  dsRed_inds <- gene_dsRed == "LentiRG"
  nondsRed_inds <- gene_dsRed == "other"
  cell_dsRedstat <- tibble(n_dsRed_umi = Matrix::colSums(seurat.list[[i]][dsRed_inds,]),
                           n_other_umi = Matrix::colSums(seurat.list[[i]][nondsRed_inds,]),
                           tot_umi = Matrix::colSums(seurat.list[[i]]),
                           prop_dsRed = n_dsRed_umi / tot_umi,
                           prop_other = n_other_umi / tot_umi,
                           SampleType = seurat.list[[i]]$SampleType,
                           Count = seurat.list[[i]]$Count,
                           Fluorescence = seurat.list[[i]]$Fluorescence)
  
  cell_dsRedstat <- cell_dsRedstat %>% mutate(dsRed.expr = case_when(n_dsRed_umi >= 1 ~ "LentiRG",TRUE ~ "No expression"))
  table3<-cell_dsRedstat %>% dplyr::count(dsRed.expr) %>% mutate(proportion = n / ncol(seurat.list[[i]]))
  show(ggplot(cell_dsRedstat, aes(n_dsRed_umi, n_other_umi, color = dsRed.expr,shape=Fluorescence)) +
         geom_point()+scale_x_log10()+scale_y_log10()+labs(title = paste0(unique(cell_dsRedstat$Count),"- LentiRG UMIs"))+scale_color_manual(values = dsRed.cols)+theme_bw())
  
  seurat.list[[i]] <- AddMetaData(seurat.list[[i]], metadata = cell_dsRedstat$dsRed.expr, col.name = "LentiRG.expression")
  tbl3<-tableGrob(table3,rows=NULL,theme=tt)
  grid.arrange(ggplot(seurat.list[[i]]@meta.data,aes(x=Count, fill=LentiRG.expression)) + geom_bar(position="fill")+
                 scale_y_continuous(labels = scales::percent)+theme_bw()+theme(axis.text.x = element_text(hjust=1,angle=45)), tbl3,ncol=2,as.table=TRUE,heights=c(2,1))
}

#QC - Part 1####
#We base the 1st QC round on nFeature_RNA count and mt% thresholds
datasetsUppertbl<-data.frame(datasetsUpper)
datasetsUppertbl$Count<-rownames(datasetsUppertbl)
datasetsHumanMitotbl<-data.frame(datasetsHumanMito)
datasetsHumanMitotbl$Count<-rownames(datasetsHumanMitotbl)
datasetsUppertbl<-inner_join(datasetsUppertbl,datasetsHumanMitotbl,by="Count")

for(i in 1:length(seurat.list)) {
  seurat.info <- tibble(Count = seurat.list[[i]]$Count)
  seurat.info <- inner_join(seurat.info,datasetsUppertbl,by="Count") 
  seurat.list[[i]]$mitoHuCH38Percent <- PercentageFeatureSet(seurat.list[[i]], pattern='*-MT-') #Calculating human mitochondrial gene percentage
  seurat.list[[i]]$mitoMM10Percent <- PercentageFeatureSet(seurat.list[[i]], pattern='*-mt-') #Calculating mouse mitochondrial gene percentage
  seurat.list[[i]] <- subset(seurat.list[[i]],subset= nFeature_RNA > 200 & nFeature_RNA < unique(seurat.info$datasetsUpper))
  seurat.list[[i]] <- subset(seurat.list[[i]],subset= mitoHuCH38Percent < unique(seurat.info$datasetsHumanMito))
  p1<-VlnPlot(seurat.list[[i]], features=c("nCount_RNA","mitoMM10Percent","mitoHuCH38Percent","nFeature_RNA"),ncol=4,group.by = "Count")+labs(title = paste0(unique(seurat.list[[i]]$Count),"-after-1st.QC"))
  p2<-FeatureScatter(seurat.list[[i]], feature1 = "nCount_RNA", feature2 = "mitoHuCH38Percent",group.by = "Count")+labs(title = paste0(unique(seurat.list[[i]]$Count),"-after-1st.QC"))
  p3<-FeatureScatter(seurat.list[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "Count")+labs(title = paste0(unique(seurat.list[[i]]$Count),"-after-1st.QC"))
  show(ggarrange(p1,p2,p3,ncol=2,nrow = 2))
}

#Identifying Doublets - DoubletFinder####

#This workflow follows the basic workflow for this package as described on their github
#(https://github.com/chris-mcginnis-ucsf/DoubletFinder)
#This is based on no ground-thruth (Because no multiplexing was performed)
#Therefore the assumed Doublet formation rate, depends on your 10X protocol and the number of expected cells
#e.g we used 5% since the cell -recovery rate for Namita's datasets was ~50.0000 and 10X gives a ~0.8% for 1000 cells
#So here you need to put in the appropriate expected multiplet rate for your dataset

seurat.list.copy<-seurat.list

for(i in 1:length(seurat.list.copy)){
  seurat.list.copy[[i]] <- SCTransform(seurat.list.copy[[i]])
  seurat.list.copy[[i]] <- RunPCA(seurat.list.copy[[i]])
  seurat.list.copy[[i]] <- FindNeighbors(object = seurat.list.copy[[i]], dims = 1:30)
  seurat.list.copy[[i]] <- FindClusters(seurat.list.copy[[i]]) #For this we will use the default resolution settings
  seurat.list.copy[[i]] <- RunUMAP(seurat.list.copy[[i]], dims = 1:30)
  
  ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
  sweep.res.seurat.list <- paramSweep_v3(seurat.list.copy[[i]], PCs = 1:ncol(seurat.list.copy[[i]]@reductions$pca) , sct = TRUE)
  sweep.stats.seurat.list <- summarizeSweep(sweep.res.seurat.list, GT = FALSE)
  bcmvn.sweep.stats.seurat.list <- find.pK(sweep.stats.seurat.list)
  
  #The default generated plot is generally really unhelpful so we generate our own
  #show(ggplot(bcmvn.sweep.stats.seurat.list, aes(pK, BCmetric, group = 1))+geom_point()+geom_line()+labs(title = paste0(dplyr::first(seurat.list.copy[[i]]$Count)," - pKs"))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  
  #This line whill then automatically select optimal pK with the hightes BCmetric
  pK <- bcmvn.sweep.stats.seurat.list %>%
    filter(BCmetric == max(BCmetric)) %>%
    select(pK) 
  pK <- as.numeric(as.character(pK[[1]]))
  
  
  
  ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
  nExp_poi <- round(0.05*nrow(seurat.list.copy[[i]]@meta.data))  ## Assuming 5% doublet formation rate - tailor for your dataset (Namita's datasets are) - This rate is technically individual for each datasets depending on the expected recovery rate during the 10X protocol
  
  annotations <- seurat.list.copy[[i]]@meta.data$seurat_clusters #based on the clustering results we're 
  homotypic.prop <- modelHomotypic(annotations)           
  
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
  seurat.list.copy[[i]] <- doubletFinder_v3(seurat.list.copy[[i]],
                                            PCs = 1:30,
                                            pN = 0.25,
                                            pK = pK,
                                            nExp = nExp_poi.adj, 
                                            reuse.pANN = FALSE,
                                            sct = TRUE)
  rm(sweep.res.seurat.list,sweep.stats.seurat.list,bcmvn.sweep.stats.seurat.list)
  
  #Showing Doublets
  p1<-DimPlot(seurat.list.copy[[i]],group.by = dplyr::first(str_subset(colnames(seurat.list.copy[[i]]@meta.data),pattern = "DF.classifications")),label=TRUE,repel=10,label.size=3)+NoLegend()+labs(title = paste0(dplyr::first(str_subset(colnames(seurat.list.copy[[i]]@meta.data),pattern = "DF.classifications"))))
  p2<-DimPlot(seurat.list.copy[[i]],group.by = dplyr::last(str_subset(colnames(seurat.list.copy[[i]]@meta.data),pattern = "DF.classifications")),label=TRUE,repel=10,label.size=3)+NoLegend()+labs(title = paste0(dplyr::first(seurat.list.copy[[i]]$Count)))
  table(pbmc.seurat.filtered@meta.data$dplyr::last(str_subset(colnames(seurat.list.copy[[i]]@meta.data),pattern = "DF.classifications")))
  p3<-DimPlot(seurat.list.copy[[i]],group.by = species,label=TRUE,repel=10,label.size=3)+NoLegend()+labs(title = "Human and Mouse")
  p4<-DimPlot(seurat.list.copy[[i]],label=TRUE,repel=10,label.size=3)+NoLegend()+labs(title = "Clusters")
  
  show(ggarrange(p1,p2,p3,p4,ncol=2,nrow=2))
  
  ##How many Doublets were found?
  print(paste0(dplyr::first(seurat.list.copy[[i]]$Count)," had ",nExp_poi.adj," doublets"))
  
  #Rename the column to get uniform meta.data columns
  
  Indx<-str_which(colnames(seurat.list.copy[[i]]@meta.data),pattern="DF.classifications")
  
  colnames(seurat.list.copy[[i]]@meta.data)[Indx] <- "Doublet/Singlet"
}

rm(seurat.list.copy)

#QC - Part 2#####
#Now we're removing those cells with species "mouse" and "mixed"

for(i in 1:length(seurat.list)){
  seurat.list[[i]] <- subset(seurat.list[[i]],subset = species == "human")
  p1<-VlnPlot(seurat.list[[i]], features=c("nCount_RNA","mitoMM10Percent","mitoHuCH38Percent","nFeature_RNA"),ncol=4,group.by = "Count")+labs(title = paste0(unique(seurat.list[[i]]$Count),"-after-2nd.QC"))
  p2<-FeatureScatter(seurat.list[[i]], feature1 = "nCount_RNA", feature2 = "mitoHuCH38Percent",group.by = "Count")+labs(title = paste0(unique(seurat.list[[i]]$Count),"-after-2nd.QC"))
  p3<-FeatureScatter(seurat.list[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "Count")+labs(title = paste0(unique(seurat.list[[i]]$Count),"-after-2nd.QC"))
  show(ggarrange(p1,p2,p3,ncol=2,nrow = 2))
  table<-seurat.list[[i]]@meta.data%>%dplyr::count(species)
  tbl<-tableGrob(table,rows=NULL,theme=tt)
  grid.arrange(ggplot(seurat.list[[i]]@meta.data,aes(x=Count, fill=species)) + geom_bar(position="fill")+
                 scale_y_continuous(labels = scales::percent)+theme_bw()+theme(axis.text.x = element_text(hjust=1,angle=45)), tbl,ncol=2,as.table=TRUE,heights=c(2,1))
}

#Preparing for multimodal reference mapping - Removing the "GRCh38-" gene prefix#####

for(i in 1:length(seurat.list)){
  updated.gene.names <- str_remove_all(rownames(seurat.list[[i]]@assays$RNA),"GRCh38-")
  seurat.list[[i]] <- RenameGenesseurat(obj = seurat.list[[i]], newnames = updated.gene.names)
}

#Multimodal reference mapping####
#In thispart, we map our datasets to a CITE-seq reference of 162,000 PBMC measured with 228 antibodies released by Satija Lab. 
#https://satijalab.org/seurat/articles/multimodal_reference_mapping.html (the vignette explaining the steps and the link to the .h5 file)
#https://www.sciencedirect.com/science/article/pii/S0092867421005833?via%3Dihub (the link to the paper)
#
#
#Merging based on SampleType and mm-reference mapping#####

#Make a safety copy of the untransformed unmerged list of seurat.objects
seurat.list.untrans.individual <- seurat.list

for(i in 1:length(seurat.list)){
  assign(unique(seurat.list[[i]]$Count),seurat.test <- seurat.list[[i]])
}

batch1list<-ls(pattern=datasetsBatch1)
batch1<-lapply(batch1list,get) #15dpi counts #n=7
batch1<-Merge_Seurat_List(batch1,add.cell.ids = NULL)

batch2list<-ls(pattern=datasetsBatch2)
batch2<-lapply(batch2list,get) #10dpt+29dpt n=4
batch2<-Merge_Seurat_List(batch2,add.cell.ids = NULL)

batch3list<-ls(pattern=datasetsBatch3)
batch3<-lapply(batch3list,get) #Uninfected n=4
batch3<-Merge_Seurat_List(batch3,add.cell.ids = NULL)

rm(batch1list,batch2list,batch3list)

list <-ls(pattern="^batch")

seurat.list<-lapply(list,get)

remove(list=ls(pattern = "^batch"))
remove(list=ls(pattern = "_count"))

###SCT Transform each Merge separately
seurat.list <- lapply(X = seurat.list, FUN = SCTransform)

#Loading the citeseq - reference
reference <- LoadH5Seurat(file = referencePath)

#Find transfer anchors between reference dataset and our seurat object
#Mapping the cell identities from the reference to the seurat object

for(i in 1:length(seurat.list)){
  if(ncol(seurat.list[[i]]) < 50){print("Less than 50 cells - not enough cells for multimodal reference mapping")}
  else{
    anchors <- FindTransferAnchors(reference = reference,query = seurat.list[[i]],normalization.method = "SCT",reference.reduction = "spca",query.assay = "RNA",dims = 1:50)
    seurat.list[[i]] <- MapQuery(anchorset = anchors,query = seurat.list[[i]],reference = reference,refdata = list(celltype.l1 = "celltype.l1",celltype.l2 = "celltype.l2",predicted_ADT = "ADT"),reference.reduction = "spca",reduction.model = "wnn.umap")
    rm(anchors)
    show(DimPlot(seurat.list[[i]], reduction = "ref.umap", group.by = "predicted.celltype.l1", label = TRUE, label.size = 3, repel = TRUE, shape.by = "Count"))
    show(DimPlot(seurat.list[[i]], reduction = "ref.umap", group.by = "predicted.celltype.l2", label = TRUE, label.size = 3, repel = TRUE, shape.by = "Count"))
    table1<-seurat.list[[i]]@meta.data%>%dplyr::count(predicted.celltype.l1)
    table2<-seurat.list[[i]]@meta.data%>%dplyr::count(predicted.celltype.l2)
    tbl1<-tableGrob(table1,rows=NULL,theme=tt)
    tbl2<-tableGrob(table2,rows=NULL,theme=tt)
    grid.arrange(ggplot(seurat.list[[i]]@meta.data,aes(x=Count, fill=predicted.celltype.l1)) + geom_bar(position="fill")+
                   scale_y_continuous(labels = scales::percent)+theme_bw()+theme(axis.text.x = element_text(hjust=1,angle=45)), tbl1,ncol=2,as.table=TRUE,heights=c(2,1))
    grid.arrange(ggplot(seurat.list[[i]]@meta.data,aes(x=Count, fill=predicted.celltype.l2)) + geom_bar(position="fill")+
                   scale_y_continuous(labels = scales::percent)+theme_bw()+theme(axis.text.x = element_text(hjust=1,angle=45)), tbl2,ncol=2,as.table=TRUE,heights=c(2,1))
  }
}


rm(reference) #Its a big object and we don't need it anymore beyond this point

#Extract predicted.celltypes and add to untr. objects####

#We're now merging all batches int one big SeuratObject to extract all the predicted.celltypes at the same time
#Seurat will automatically add individual suffixes to the cells which we will have to remove later on

seuratrefmapped.merged<-Merge_Seurat_List(list_seurat = seurat.list,add.cell.ids = NULL)

celltypes.after.mm <- FetchData(seuratrefmapped.merged, vars = c('predicted.celltype.l1', 'predicted.celltype.l2','Count','SampleType','predicted.celltype.l1.score','predicted.celltype.l2.score'), slot="counts")

rm(seuratrefmapped.merged)

#We split the celltypes.after.mm table and set the CellIDs as rownames -> we can't use the default rownames, since those are with added suffixes
#so we need to remove that first before we can add metadata to the untransformed datasets

celltypes.after.mm.split <-split(celltypes.after.mm,celltypes.after.mm$Count)
for (i in 1:length(celltypes.after.mm.split)) {
  celltypes.after.mm.split[[i]]$cell.id<-substr(rownames(celltypes.after.mm.split[[i]]),1,18)
  rownames(celltypes.after.mm.split[[i]]) <- celltypes.after.mm.split[[i]]$cell.id
}

##Adds the predicted.celltype.l1 and l2 as metadata (technically we also add Count and SampleType, but those just get replaced with the same thing)

for (i in 1:length(seurat.list.untrans.individual)) {
  seurat.list.untrans.individual[[i]]<-AddMetaData(seurat.list.untrans.individual[[i]],celltypes.after.mm.split[[i]])
}


##The non-normalized individual datasets now include predicted celltypes as meta.data
##Based on mapping to the cite-seq reference from the Satija Lab reference dataset
##For the rest of the workflow we chose to merge datasets into seurat.list based on mouse/or type of sample type
#thereby minimizing batch effects

#Normalization####
#SCT transform the batches

for (i in 1:length(seurat.list)) {
  seurat.list[[i]]<-SCTransform(seurat.list[[i]],vars.to.regress = c("mitoHuCH38Percent","mitoMM10Percent"))
}

#JackStraw and Elbow Plot

#JackStraw and 1st.PCA (needed to decide on thresholds)####
#If we want to do the JackStraw and all that
#JackStraw
seurat.list.copy<-seurat.list


for (i in 1:length(seurat.list.copy)){
    DefaultAssay(seurat.list.copy[[i]]) <- "RNA"
    print("Performing JackStraw and PCA")
    #Normalize RNA-assay data
    seurat.list.copy[[i]]<-NormalizeData(seurat.list.copy[[i]])
    
    #Scale RNA assay data
    seurat.list.copy[[i]]<-ScaleData(seurat.list.copy[[i]])
    
    #Find variabke features
    seurat.list.copy[[i]]<-FindVariableFeatures(seurat.list.copy[[i]])
    
    test.features <- VariableFeatures(seurat.list.copy[[i]])
    test.features <- setdiff(test.features,featuresHIV)
    mouse.genes<-str_subset(test.features,pattern = "mm10")
    test.features <- setdiff(test.features,mouse.genes)
    test.features <- setdiff(test.features,featuresUCOERG)
    
    #RUN PCA
    if (ncol(seurat.list.copy[[i]]@assays$RNA)<500){
      print("Datasest contains less than 500 cells - RunPCA with npcs=15")
      seurat.list.copy[[i]]<-RunPCA(seurat.list.copy[[i]],features=test.features,npcs = 15)
      print("Datasest contains less than 500 cells - JackStraw with dim=15")
      seurat.list.copy[[i]] <- JackStraw(seurat.list.copy[[i]], num.replicate =100,dims=15)
      seurat.list.copy[[i]] <- ScoreJackStraw(seurat.list.copy[[i]], dims =1:15)
      show(JackStrawPlot(seurat.list.copy[[i]], dims= 1:15)+
             labs(title = paste0(dplyr::first(seurat.list.copy[[i]]$Count)," n=",ncol(seurat.list.copy[[i]]@reductions$pca))))
      show(ElbowPlot(seurat.list.copy[[i]],ndims=15)+
             labs(title = paste0(dplyr::first(seurat.list.copy[[i]]$Count)," n=",ncol(seurat.list.copy[[i]]@reductions$pca))))
      
    }
    
    else if (ncol(seurat.list.copy[[i]]@assays$RNA) %in% c(500:5000)){
      print("Dataset contains more than 500 but less than 5000 cells - RunPCA with npcs=35")
      seurat.list.copy[[i]]<-RunPCA(seurat.list.copy[[i]],features=test.features,npcs = 35)
      print("Dataset contains more than 500 but less than 5000 cells - JackStraw with dims=30")
      seurat.list.copy[[i]] <- JackStraw(seurat.list.copy[[i]], num.replicate =100,dims=30)
      seurat.list.copy[[i]] <- ScoreJackStraw(seurat.list.copy[[i]], dims =1:30)
      show(JackStrawPlot(seurat.list.copy[[i]], dims= 1:30)+
             labs(title = paste0(dplyr::first(seurat.list.copy[[i]]$Count)," n=",ncol(seurat.list.copy[[i]]@reductions$pca))))
      show(ElbowPlot(seurat.list.copy[[i]],ndims=30)+
             labs(title = paste0(dplyr::first(seurat.list.copy[[i]]$Count)," n=",ncol(seurat.list.copy[[i]]@reductions$pca))))
    }
    else if (ncol(seurat.list.copy[[i]]@assays$RNA)>5000){
      print("Dataset contains more than 5000 cells - RunPCA with npcs=45")
      seurat.list.copy[[i]]<-RunPCA(seurat.list.copy[[i]],features=test.features,npcs = 45)
      print("Dataset contains more than 5000 cells - JackStraw with dims=45")
      seurat.list.copy[[i]] <- JackStraw(seurat.list.copy[[i]], num.replicate =100,dims=45)
      seurat.list.copy[[i]] <- ScoreJackStraw(seurat.list.copy[[i]], dims =1:45)
      show(JackStrawPlot(seurat.list.copy[[i]], dims= 1:45)+
             labs(title = paste0(dplyr::first(seurat.list.copy[[i]]$Count)," n=",ncol(seurat.list.copy[[i]]@reductions$pca))))
      show(ElbowPlot(seurat.list.copy[[i]],ndims=45)+
             labs(title = paste0(dplyr::first(seurat.list.copy[[i]]$Count)," n=",ncol(seurat.list.copy[[i]]@reductions$pca))))
    } 
    
}

#The chosen thresholds can be found in the supplementary material


#Clustree - Determine cluster resolution####
#Now we're doing the PCA/UMAP/TSNE based on the number of dimensions we determined using the JackStraw/ElbowPlot
#And look at the clustree plot to determine the correct resolution for the UMAPs

#https://lazappi.github.io/clustree/articles/clustree.html
#This link will show the documentation for the resulting plot
#Based on this you can decide the optimal resolution when clustering, to not over/under cluster your cells

datasetsUMAPtbl<-data.frame(datasetsUMAP)
datasetsUMAPtbl$Count<-rownames(datasetsUMAPtbl)

resolution.range <- seq(from = 0, to = 1, by = 0.1)

seurat.list.copy<-seurat.list

for (i in 1:length(seurat.list.copy)){
  seurat.info <- tibble(Count = seurat.list.copy[[i]]$Count)
  seurat.info <- inner_join(seurat.info,datasetsUMAPtbl,by="Count")
    
  test.features <- VariableFeatures(seurat.list.copy[[i]])
  test.features <- setdiff(test.features,featuresHIV)
  mouse.genes<-str_subset(test.features,pattern = "mm10")
  test.features <- setdiff(test.features,mouse.genes)
  test.features <- setdiff(test.features,featuresUCOERG)
  seurat.list.copy[[i]] <- RunPCA(seurat.list.copy[[i]],features=test.features)
  seurat.list.copy[[i]] <- RunUMAP(seurat.list.copy[[i]],dims = 1:unique(seurat.info$datasetsUMAP))
  seurat.list.copy[[i]] <- RunTSNE(seurat.list.copy[[i]],dims = 1:unique(seurat.info$datasetsUMAP))
  seurat.list.copy[[i]] <- FindNeighbors(seurat.list.copy[[i]],dims = 1:unique(seurat.info$datasetsUMAP))
  seurat.list.copy[[i]] <- FindClusters(seurat.list.copy[[i]],resolution = resolution.range) 
    
  show(clustree(seurat.list.copy[[i]], prefix = "SCT_snn_res."))
}


#Based on the chosen resolution, which can be found in the supplementary data,
#we create the UMAP and Heatmap plots (pre-subsetting)

datasetsClusterRestbl<-data.frame(datasetsClusterRes)
datasetsClusterRestbl$Count<-rownames(datasetsClusterRestbl)

datasetsUMAPandClusterRestbl<-inner_join(datasetsUMAPtbl,datasetsClusterRestbl,by="Count")


for (i in 1:length(seurat.list)){
    
    seurat.info <- tibble(Count = seurat.list[[i]]$Count)
    seurat.info <- inner_join(seurat.info,datasetsUMAPandClusterRestbl,by="Count")
    
    test.features <- VariableFeatures(seurat.list[[i]])
    test.features <- setdiff(test.features,featuresHIV)
    mouse.genes<-str_subset(test.features,pattern = "mm10")
    test.features <- setdiff(test.features,mouse.genes)
    test.features <- setdiff(test.features,featuresUCOERG)
    seurat.list[[i]] <- RunPCA(seurat.list[[i]],features=test.features)
    seurat.list[[i]] <- RunUMAP(seurat.list[[i]],dims = 1:unique(seurat.info$datasetsUMAP))
    seurat.list[[i]] <- RunTSNE(seurat.list[[i]],dims = 1:unique(seurat.info$datasetsUMAP))
    seurat.list[[i]] <- FindNeighbors(seurat.list[[i]],dims = 1:unique(seurat.info$datasetsUMAP))
    seurat.list[[i]] <- FindClusters(seurat.list[[i]],resolution = unique(seurat.info$datasetsClusterRes)) 
    
    p1<-DimPlot(seurat.list[[i]],label = TRUE,repel = 10,label.size = 3)+NoLegend()+labs(title = paste0(dplyr::first(seurat.list[[i]]$Count),"-UMAP - pre-subsetting"),subtitle = paste0("ndims=1:",dplyr::first(seurat.info$datasetsUMAP)," res=",dplyr::first(seurat.info$datasetsClusterRes)))
    p2<-DimPlot(seurat.list[[i]],group.by = "Count",label = TRUE,repel = 10,label.size = 3)+NoLegend()+labs(title = paste0(dplyr::first(seurat.list[[i]]$Count),"-UMAP - pre-subsetting"),subtitle = paste0("ndims=1:",dplyr::first(seurat.info$datasetsUMAP)," res=",dplyr::first(seurat.info$datasetsClusterRes)))
    p3<-DimPlot(seurat.list[[i]],group.by = "predicted.celltype.l1",label = TRUE,repel = 10,label.size = 3)+NoLegend()+labs(title = paste0(dplyr::first(seurat.list[[i]]$Count),"-UMAP - pre-subsetting"),subtitle = paste0("ndims=1:",dplyr::first(seurat.info$datasetsUMAP)," res=",dplyr::first(seurat.info$datasetsClusterRes)))
    p4<-DimPlot(seurat.list[[i]],group.by = "predicted.celltype.l2",label = TRUE,repel = 10,label.size = 3)+NoLegend()+labs(title = paste0(dplyr::first(seurat.list[[i]]$Count),"-UMAP - pre-subsetting"),subtitle = paste0("ndims=1:",dplyr::first(seurat.info$datasetsUMAP)," res=",dplyr::first(seurat.info$datasetsClusterRes)))
    show(ggarrange(p1,p2,p3,p4,ncol=2,nrow = 2))
    
    p5<-DimPlot(seurat.list[[i]],group.by = "Doublet.Singlet",label = TRUE,repel = 10,label.size = 3)+NoLegend()+labs(title = paste0(dplyr::first(seurat.list[[i]]$Count),"-UMAP - pre-subsetting"),subtitle = paste0("ndims=1:",dplyr::first(seurat.info$datasetsUMAP)," res=",dplyr::first(seurat.info$datasetsClusterRes)))
    p6<-DimPlot(seurat.list[[i]],group.by = "status",label = TRUE,repel = 10,label.size = 3,order=c("HIV+ very high","HIV+ high","HIV+ low","HIV-"),cols = HIV.cols)+NoLegend()+labs(title = paste0(dplyr::first(seurat.list[[i]]$Count),"-UMAP - pre-subsetting"),subtitle = paste0("ndims=1:",dplyr::first(seurat.info$datasetsUMAP)," res=",dplyr::first(seurat.info$datasetsClusterRes)))
    p7<-DimPlot(seurat.list[[i]],group.by = "Fluorescence",label = TRUE,repel = 10,label.size = 3,order=c("GFP","dsRed","Mixed","Unmarked"),cols = fluor.cols)+NoLegend()+labs(title = paste0(dplyr::first(seurat.list[[i]]$Count),"-UMAP - pre-subsetting"),subtitle = paste0("ndims=1:",dplyr::first(seurat.info$datasetsUMAP)," res=",dplyr::first(seurat.info$datasetsClusterRes)))
    p8<-DimPlot(seurat.list[[i]],group.by = "LentiRG.expression",label = TRUE,repel = 10,label.size = 3,order=c("LentiRG","No epxression"),cols = dsRed.cols)+NoLegend()+labs(title = paste0(dplyr::first(seurat.list[[i]]$Count),"-UMAP - pre-subsetting"),subtitle = paste0("ndims=1:",dplyr::first(seurat.info$datasetsUMAP)," res=",dplyr::first(seurat.info$datasetsClusterRes)))
    show(ggarrange(p6,p7,p8,ncol=2,nrow = 2))
    
    ##Barplots
    
    a1<-ggplot(seurat.list[[i]]@meta.data,aes(x=seurat_clusters, fill=predicted.celltype.l1)) + geom_bar(position="fill")+
      scale_y_continuous(labels = scales::percent)+theme_bw()+labs(title = paste0(dplyr::first(seurat.list[[i]]$Count)),
                                                                   subtitle = "Before cluster subsetting")+
      theme(axis.text.x = element_text(hjust=1,angle=45),legend.key.size = unit(0.1, 'cm'))
    
    a2<-ggplot(seurat.list[[i]]@meta.data,aes(x=seurat_clusters, fill=predicted.celltype.l2)) + geom_bar(position="fill")+
      scale_y_continuous(labels = scales::percent)+theme_bw()+labs(title = paste0(dplyr::first(seurat.list[[i]]$Count)),
                                                                   subtitle = "Before cluster subsetting")+
      theme(axis.text.x = element_text(hjust=1,angle=45),legend.key.size = unit(0.1, 'cm'))
    show(ggarrange(a1,a2,ncol=1,nrow = 2))
    
    h1<-ggplot(seurat.list[[i]]@meta.data,aes(x=seurat_clusters, fill=status)) + geom_bar(position="fill")+
      scale_y_continuous(labels = scales::percent)+theme_bw()+labs(title = paste0(dplyr::first(seurat.list[[i]]$Count)),
                                                                   subtitle = "Before cluster subsetting")+
      theme(axis.text.x = element_text(hjust=1,angle=45),legend.key.size = unit(0.1, 'cm'))+
      scale_fill_manual(values = HIV.cols)
    
    h2<-ggplot(seurat.list[[i]]@meta.data,aes(x=predicted.celltype.l1, fill=status)) + geom_bar(position="fill")+
      scale_y_continuous(labels = scales::percent)+theme_bw()+labs(title = paste0(dplyr::first(seurat.list[[i]]$Count)),
                                                                   subtitle = "Before cluster subsetting")+
      theme(axis.text.x = element_text(hjust=1,angle=45),legend.key.size = unit(0.1, 'cm'))+
      scale_fill_manual(values = HIV.cols)
    
    show(ggarrange(h1,h2,ncol=1,nrow = 2))
    
    tableXYZ<-seurat.list[[i]]@meta.data%>%
      group_by(seurat_clusters)%>%
      dplyr::count(status)%>%
      mutate(percent=(n/sum(n))*100)
    
    tt<-ttheme_default(colhead=list(fg_params=list(parse=TRUE)),base_size = 5)
    tblxyz<-tableGrob(tableXYZ,rows=NULL,theme=tt)
    
    show(grid.arrange(h1, tblxyz,
                      ncol=2,
                      as.table=TRUE))
    
    t1<-ggplot(seurat.list[[i]]@meta.data,aes(x=seurat_clusters, fill=LentiRG.expression)) + geom_bar(position="fill")+
      scale_y_continuous(labels = scales::percent)+theme_bw()+labs(title = paste0(dplyr::first(seurat.list[[i]]$Count)),
                                                                   subtitle = "Before cluster subsetting")+
      theme(axis.text.x = element_text(hjust=1,angle=45),legend.key.size = unit(0.1, 'cm'))+
      scale_fill_manual(values = dsRed.cols)
    
    t2<-ggplot(seurat.list[[i]]@meta.data,aes(x=predicted.celltype.l1, fill=LentiRG.expression)) + geom_bar(position="fill")+
      scale_y_continuous(labels = scales::percent)+theme_bw()+labs(title = paste0(dplyr:first(seurat.list[[i]]$Count)),
                                                                   subtitle = "Before cluster subsetting")+
      theme(axis.text.x = element_text(hjust=1,angle=45),legend.key.size = unit(0.1, 'cm'))+
      scale_fill_manual(values = dsRed.cols)
    
    show(ggarrange(t1,t2,ncol=1,nrow = 2))
    
    tableXYZ<-seurat.list[[i]]@meta.data%>%
      group_by(seurat_clusters)%>%
      dplyr::count(LentiRG.expression)%>%
      mutate(percent=(n/sum(n))*100)
    
    tblxyz<-tableGrob(tableXYZ,rows=NULL,theme=tt)
    
    show(grid.arrange(t1, tblxyz,
                      ncol=2,
                      as.table=TRUE))
    
    o1<-ggplot(seurat.list[[i]]@meta.data,aes(x=seurat_clusters, fill=Fluorescence)) + geom_bar(position="fill")+
      scale_y_continuous(labels = scales::percent)+theme_bw()+labs(title = paste0(dplyr::first(seurat.list[[i]]$Count)),
                                                                   subtitle = "Before cluster subsetting")+
      theme(axis.text.x = element_text(hjust=1,angle=45),legend.key.size = unit(0.1, 'cm'))+
      scale_fill_manual(values = fluor.cols)
    
    o2<-ggplot(seurat.list[[i]]@meta.data,aes(x=predicted.celltype.l1, fill=Fluorescence)) + geom_bar(position="fill")+
      scale_y_continuous(labels = scales::percent)+theme_bw()+labs(title = paste0(dplyr::first(seurat.list[[i]]$Count)),
                                                                   subtitle = "Before cluster subsetting")+
      theme(axis.text.x = element_text(hjust=1,angle=45),legend.key.size = unit(0.1, 'cm'))+
      scale_fill_manual(values = fluor.cols)
    
    show(ggarrange(o1,o2,ncol=1,nrow = 2))
    
    tableXYZ<-seurat.list[[i]]@meta.data%>%
      group_by(seurat_clusters)%>%
      dplyr::count(Fluorescence)%>%
      mutate(percent=(n/sum(n))*100)
    
    tblxyz<-tableGrob(tableXYZ,rows=NULL,theme=tt)
    
    show(grid.arrange(o1, tblxyz,
                      ncol=2,
                      as.table=TRUE))
    
    #Heatmaps presubsetting
    seurat.list.markers <- FindAllMarkers(seurat.list[[i]])
    mousemarkers<-str_subset(seurat.list.markers$gene,pattern = "mm10")
    seurat.list.markers<-seurat.list.markers[!seurat.list.markers$gene %in% mousemarkers,]
    seurat.list.markers<-seurat.list.markers[!seurat.list.markers$gene %in% featuresHIV,]
    seurat.list.markers<-seurat.list.markers[!seurat.list.markers$gene %in% featuresUCOERG,]
    
    #Save Markers full list and Top10
    assign(paste0(dplyr::first(seurat.list[[i]]$Count),"-markers-unsupervised"),seurat.list.markers)
    seurat.list.markers.top10 <- seurat.list.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
    assign(paste0(dplyr::first(seurat.list[[i]]$Count),"-markers-unsupervised-topten"),seurat.list.markers.top10)
    
    #Show Heatmap
    show(DoHeatmap(seurat.list[[i]], features = seurat.list.markers.top10$gene,size = 3)+NoLegend()+
           labs(title = paste0(dplyr::first(seurat.list[[i]]$Count)),
                subtitle = "Before cluster subsetting")+theme(text = element_text(size = 5)))
    
}

#Finish pdf#####
dev.off()

#Save workspace and seurat.list and move on to Workflow-canonical-cellmarkers.R

save.image("Seurat.RData")
saveRDS(seurat.list,"seurat.list.pre.subset.rds")





