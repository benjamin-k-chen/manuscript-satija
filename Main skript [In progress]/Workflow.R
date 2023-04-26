#Seurat downstream analysis - BASIC WORKFLOW - 2023-03
#
#
#This is the basic workflow for the analysis of the paper Satija et al. 2023
#
#This workflow was performed on each data set separately until MMRM
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

#To be able to accurately do perform the multimodal reference mapping the gene names in both the reference and our dataset should match
#Since we used the 10X cellranger reference with both mouse and human genome we have to take the prefix out of the dataset
#We will create the function "RenameGeneseurat" that replaces the gene names (rownames within defined slots of the object)

RenameGenesseurat <- function(obj = seurat.object, newnames = updated.gene.names) { # Replace gene names in different slots of a seurat object. Run this before any type of integration. It only changes obj@assays$RNA@counts, @data and @scale.data.
  print("Run this before integration. It only changes obj@assays$RNA@counts, @data")
  RNA <- obj@assays$RNA
  if (nrow(RNA) == length(newnames)) {
    if (length(RNA@counts)) RNA@counts@Dimnames[[1]]            <- newnames
    if (length(RNA@data)) RNA@data@Dimnames[[1]]                <- newnames
  } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
  obj@assays$RNA <- RNA
  return(obj)
}

#Apply the function to our SeuratObject
updated.gene.names <- str_remove_all(rownames(seurat.object@assays$RNA),"GRCh38-") #Making a vector of the gene names without the prefix

#Applying the function
seurat.object <- RenameGenesseurat(obj = seurat.object, newnames = updated.gene.names)

#Multimodal reference mapping####
#In thispart, we map our datasets to a CITE-seq reference of 162,000 PBMC measured with 228 antibodies released by Satija Lab. 
#https://satijalab.org/seurat/articles/multimodal_reference_mapping.html (the vignette explaining the steps and the link to the .h5 file)
#https://www.sciencedirect.com/science/article/pii/S0092867421005833?via%3Dihub (the link to the paper)
#
#
#Similar to DoubletFinder we will first create a copy of our object to which we will later add the predicted.celltypes
seurat.object.original<-seurat.object

#SCT normalizing the dataset
seurat.object <- SCTransform(seurat.object)

#Loading the citeseq - reference
reference <- LoadH5Seurat(file = "/Path/to/the/h5/file")

#Find transfer anchors between reference dataset and our seurat object

anchors <- FindTransferAnchors(reference = reference,query = seurat.object,
                               normalization.method = "SCT",
                               reference.reduction = "spca",
                               query.assay = "RNA",
                               dims = 1:50) #We chose 50 dimension, but you can adjust that to your dataset as you like

#Mapping the cell identities from the reference to the seurat object
seurat.object <- MapQuery(anchorset = anchors,
                          query = seurat.object,
                          reference = reference,
                          refdata = list(celltype.l1 = "celltype.l1",celltype.l2 = "celltype.l2",predicted_ADT = "ADT"),
                          reference.reduction = "spca",
                          reduction.model = "wnn.umap")

#We can display our cells in the UMAP of the reference dataset
DimPlot(seurat.object, reduction = "ref.umap",
        group.by = "predicted.celltype.l1",
        label = TRUE, label.size = 3, repel = TRUE)
    
DimPlot(seurat.object, reduction = "ref.umap",
        group.by = "predicted.celltype.l2", 
        label = TRUE, label.size = 3, repel = TRUE)

#Let's make tables looking at the number of predicted celltypes

#Level 1 - Cell Types
table1<-seurat.object@meta.data%>%dplyr::count(predicted.celltype.l1)

#Level 2 - Cell Types
table2<-seurat.object@meta.data%>%dplyr::count(predicted.celltype.l2)

#After this we can remove the reference to clean up the R enviroment
rm(reference)

#We don't need the normalization and PCA/UMAP anymore so we take the meta.data and add it to the copy of seurat.object

Celltypes.after.mm <- FetchData(seurat.object, vars = c('predicted.celltype.l1','predicted.celltype.l2'), slot="counts") #In this table the cell.barcodes are the rownames, those are matching to the copy and are used to match the meta.data to the right cell

seurat.object.original<-AddMetaData(seurat.object.original,Celltypes.after.mm) #The automatic column name will be the column name of the table

#We can now safely remove the SeuratObject we don't need by overwriting them

seurat.object<-seurat.object.original

rm(seurat.object.original)

#Normalization####
#At this point we follow the basic Seurat workflow released by Satija Lab
#https://satijalab.org/seurat/articles/pbmc3k_tutorial.html#normalizing-the-data-1
#When doing SCT-transform we don't need to perform Normalize/ScaleData anymore (unless when performing a JackStraw)

seurat.object<-SCTransform(seurat.object,vars.to.regress = c("mitoHuCH38Percent","mitoMM10Percent"))

#Identification of highly variable features####
#We calculate a subset of features that a high high cell-to-cell variation in the dataset
#(i.e, they are highly expressed in some cells, and lowly expressed in others)
#Focusing on these genes in downstream analysis helps to highlight biological signal scRNA datasets.

seurat.object<-FindVariableFeatures(seurat.object,selection.method = "vst",nfeatures = 2000)

#Take the top10 most highly variable genes in the dataset
top10<- head(VariableFeatures(seurat.object),10)

#Then we plot these genes in the VariableFeaturePlot
plot1 <- VariableFeaturePlot(seurat.object)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 | plot2


#JackStraw - Determining the dimensionality of the dataset####

#We first make a copy of the dataset
seurat.object.original<-seurat.object

#Unfortunately the JackStraw workflow doesn't work on our SCT transformed dataset, so we need to change the default assay
#back to "RNA" (untransformed/unnormalized raw data) and normalize and scale that data using default settings

DefaultAssay(seurat.object) <- "RNA"

#Normalize RNA-assay data
seurat.object<-NormalizeData(seurat.object)
    
#Scale RNA assay data
seurat.object<-ScaleData(seurat.object)

#Find variable features
seurat.object<-FindVariableFeatures(seurat.object)

#We make a variable of the highly variable features and take out the mouse genes - HIV genes - Lenti-RG genes
#That could have an impact on the PCA and UMAP

featuresUCOERG<- c("dsRed", "EGFP", "WPRE" )
featuresHIV<- c("gag-pol","pol","pol-vif-vpr-tat-rev","vpu-env","env","mirfp670nano","p2a-cre-ires-nef" ) #This of course depends on the genes defined in your cellranger reference

#
features <- VariableFeatures(seurat.object)
features <- setdiff(features,featuresHIV) #Taking out the HIV genes in the list of highly variable genes
mouse.genes<-str_subset(features,pattern = "mm10") #Making a list of the genes that are highly variable and have the mouse gene prefix
features <- setdiff(features,mouse.genes) #Taking out the mouse genes in the list of highly variable genes
features <- setdiff(features,featuresUCOERG) #Taking out the Lenti-RG genes in the list of highly variable genes
    
#RUN a basic PCA
seurat.object<-RunPCA(seurat.object,features=features,npcs = 45)

#Performing a JackStraw
seurat.object <- JackStraw(seurat.object, num.replicate=100)
seurat.object <- ScoreJackStraw(seurat.object, dims=1:45)

#Plot the JackStraw results
JackStrawPlot(seurat.object, dims=1:45)+
       labs(title = paste0(" n=",ncol(seurat.objects[[i]]@reductions$pca)))

ElbowPlot(seurat.object,ndims=45)+
       labs(title = paste0("n=",ncol(seurat.object@reductions$pca)))

#These plots allow us to determine the dimensionality of the dataset for the later downstream analysis
#Here you need to decide on thresholds and mark them down!


#Clustree - Determine cluster resolution####
#https://lazappi.github.io/clustree/articles/clustree.html
#This link will show the documentation for the resulting plot
#Based on this you can decide the optimal resolution when clustering, to not over/under cluster your cells

#We first make a copy of the dataset
seurat.object.original<-seurat.object

#We make a range of cluster resolutions
resolution.range <- seq(from = 0, to = 1, by = 0.1)

#We then again identify the highly variable genes in the dataset and take out Lenti and HIV genes
features <- VariableFeatures(seurat.object)
features <- setdiff(features,featuresHIV)
mouse.genes<-str_subset(features,pattern = "mm10")
features <- setdiff(features,mouse.genes)
features <- setdiff(features,featuresUCOERG)

#We will have to run a basic PCA and UMAP sung the thresholds we have previously defined (->JackStraw)

seurat.object <- RunPCA(seurat.object,features=test.features)
seurat.object <- RunUMAP(seurat.object,dims = 1:XX) #Insert the threshold, in place of "XX"
seurat.object <- FindNeighbors(seurat.object,dims = 1:XX)
seurat.object <- FindClusters(seurat.object,resolution = resolution.range) 

#Now lets look at the clustertree 
clustree(seurat.object, prefix = "SCT_snn_res.")

#Now you can decide on the optimal resolution -> Write it down!

#Creating the UMAPs####
#We can't use the same previous seurat.object so we will overwrite it with the copey we made before running clustree
seurat.object<-seurat.object.original

seurat.object <- RunPCA(seurat.object,features=test.features)
seurat.object <- RunUMAP(seurat.object,dims = 1:XX) #Insert the threshold, in place of "XX"
seurat.object <- FindNeighbors(seurat.object,dims = 1:XX)
seurat.object <- FindClusters(seurat.object,resolution = YY) #Insert the resolution, in place of "YY"

#Show the UMAP plot
#We're gonna save four dfferent UMAP plots and them shwo them in one figure
#In this iteration we have the labels for the colors inside the UMAP and no legend on the right
#If you want a legend and no labels just remove ",label = TRUE,repel = 10,label.size = " and "+NoLegend()" from the code

p1<-DimPlot(seurat.object,label = TRUE,repel = 10,label.size = 3)+NoLegend()
p2<-DimPlot(seurat.object,group.by = "ANY META-DATA COLUMN",label = TRUE,repel = 10,label.size = 3)+NoLegend()
p3<-DimPlot(seurat.object,group.by = "predicted.celltype.l1",label = TRUE,repel = 10,label.size = 3)+NoLegend() 
p4<-DimPlot(seurat.object,group.by = "predicted.celltype.l2",label = TRUE,repel = 10,label.size = 3)+NoLegend()

ggarrange(p1,p2,p3,p4,ncol=2,nrow = 2)

#You can also manually define colors and use them in the display of your UMAP
#e.g.:
HIV.cols <- c('HIV-'='lightgrey','HIV+ low'="#ffd8b1",'HIV+ high'= "purple",'HIV+ very high'="purple") 
#There are predefined colors in R of which you can just type in the name or you can use custom colors using the Hexcode "#000000"
#Each name matches one of the definitions on our meta.data column
#With adding "cols=___" in DimPlot() you can define your colors
DimPlot(seurat.object,group.by = "ANY META DATA COLUMN",label = TRUE,repel = 10,label.size = 3,cols = HIV.cols)+NoLegend()

#Creating HeatMaps with top10 marker genes per cluster####
#Heatmaps are created by comparing all the cluster to one another and identifying defininig marker genes for each

seurat.object.markers <- FindAllMarkers(seurat.object) #We will then again remove the HIV and Lenti-RG genes from this list of markers
mousemarkers<-str_subset(seurat.object$gene,pattern = "mm10")
seurat.object.markers<-seurat.object.markers[!seurat.object.markers$gene %in% mousemarkers,]
seurat.object.markers<-seurat.object.markers[!seurat.object.markers$gene %in% featuresHIV,]
seurat.object.markers<-seurat.object.markers[!seurat.object.markers$gene %in% featuresUCOERG,]

#Now we're only taking the top10 markers in each cluster
seurat.object.markers.top10 <- seurat.object.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

#Show Heatmap
DoHeatmap(seurat.object, features = seurat.object.markers.top10$gene,size = 3)

#Now you can look at the heatmap with the marker genes for each cluster in your dataset
#Ideally you should have a descending stair shape
#If not you might have to go back a couple steps and troubleshoot
#This Heatmap is also incredibly helpful if you want to identify your own celltypes not based on reference mapping
#In that case you have to go through the top markers for each cluster manually and use that top10 list to find distinctive gene markers based on literature

#FeaturePlots#####
##Another way to look at your UMAPs is by looking at the expression of certain cell type markers in your cells
#We made a defined list of cellmarker genes at which we wan to look at
#You can obviously adjust that list to your liking - depending on what celltypes you expect in your dataset

#T.cell("CD3D","CD3E","CD3G","TRAC") B.cell("CD19","MS4A1") NK.cells("KLRF1") MonoCD14CD16("CD14","FCGR3A")
cellmarker.features<-c("CD3D","CD3E","CD3G","TRAC","CD19","MS4A1","KLRF1","CD14","FCGR3A")

#Plot the expression in a 'FeaturePLot'
FeaturePlot(seurat.object,features = cellmarker.features)

#You can also look at the expression levels of these genes in a violin plot

VlnPlot(seurat.object,features = cellmarker.features) #Grouped by seurat defined clusters

VlnPlot(seurat.object,features = cellmarker.features, group.by = "ANY META DATA COLUMN")

#Removing clusters or low abundance celltypes#####
#At this point you might have look at your UMAPs and decide to clean up unwanted clusters and/or low abundance cell types
#Please keep in mind though: If you remove cells you have to go all the way back to Normalization on reiterate the whole workflow process

#Example: We looked at our UMAP and HeatMap and decided that clusters 3 and 4 out of the total 9 (0-9) are not needed

seurat.object<-subset(seurat.object,subset = seurat_clusters %in% c(0,1,2,5,6,7,8,9)) #Basically here we tell Seurat to take all the clusters EXCEPT 3 and 4


#Example: We looked at our UMAP and HeatMap and decided all cells under 1% total are to be removed UNLESS they are HIV+ or GFP+ (several conditions)

#First we need to create a table of all predicted celltypes in our dataset and their number

celltypes<-seurat.object@meta.data%>%
  dplyr::count(predicted.celltype.l2)%>%
  mutate(percent=(n/sum(n))*100)

celltypes.more.than.1.percent<-celltypes[celltypes$percent > 1,] #Now we're taking all the celltypes that have a higher than >1% proportion

celltypes.more.than.1.percent.list<-celltypes.more.than.1.percent$predicted.celltype.l2 #Save all the names in a list/value

#Now we subset
#Our Conditions are seperated with "|" 
#"|" means "or", but if want to select cells that match all these criteria you use "&&" which means "and"

seurat.object <- subset(seurat.object,subset = predicted.celltype.l2 %in% celltypes.more.than.1.percent.list |
                                                         status == c("HIV+ low, HIV+ high, HIV+ very high") |
                                                         Fluorescence == "GFP")


