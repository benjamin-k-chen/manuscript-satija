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

seurat.object <- Read10X(data.dir = "path/to/filtered_feature_bc_matrix/") 
seurat.object <- CreateSeuratObject(counts = seurat.object, project = "NAME")

#Mitochondrial gene percentage#####
##Our custom reference was was made using the 10x reference that contains both human and mouse genes,
#therefore we can calculate both mitochondrial gene percentages

seurat.object$mitoHuCH38Percent <- PercentageFeatureSet(seurat.object, pattern='*-MT-') #Calculating human mitochondrial gene percentage
seurat.object$mitoMM10Percent <- PercentageFeatureSet(seurat.object, pattern='*-mt-') #Calculating mouse mitochondrial gene percentage

#Violin Plots######
#This shows the violin plot for the four categories
VlnPlot(seurat.object, features=c("nCount_RNA","mitoMM10Percent","mitoHuCH38Percent","nFeature_RNA")) 

#Scatter Plots####
FeatureScatter(seurat.object, feature1 = "nCount_RNA", feature2 = "mitoHuCH38Percent")
FeatureScatter(seurat.object, feature1 = "nCount_RNA", feature2 = "mitoMM10Percent")
FeatureScatter(seurat.object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

#The individual threshold values can be found

#Metadata######
##Add some meta.data to our dataset
#Based on what we know i.e. "SampleType","Treatment","Fluorescence","Condition" etc.
#You add as many new columns as you like

seurat.object <- AddMetaData(seurat.object, metadata = "______", col.name = "_____")


#Species Identification#### 
#If you used a the 10X reference that contains both mouse and human genes

#We define every gene that has the prefix "mm10" as "mouse" and every other as "human"
#This basically creates a vector based on the genes in the dataset (rownames) and gives them the indicator "human" or "mouse"
gene_species <- ifelse(str_detect(rownames(seurat.object),"mm10*"), "mouse", "human")

#We then use these indicators to get a True/False statement for each gene in that list wether they are "human" or "mouse"
mouse_inds <- gene_species == "mouse" #All "mouse" genes get indicator "TRUE" all "human" get "FALSE"
human_inds <- gene_species == "human" #All "human" genes get indicator "TRUE" all "mouse" get "FALSE"

#Using these indicators we can now create a table calles "cell_species" using tibble
#We calculate the amount of "human" or "mouse" umi based on the amount of genes that are "TRUE" within the total umi count
#of each cell which we calculate using colSums() 
#Within this table we also calculate the % of human/mouse umi within each cell (each row represents a cell)

cell_species <- tibble(n_mouse_umi = Matrix::colSums(seurat.object[mouse_inds,]),
                       n_human_umi = Matrix::colSums(seurat.object[human_inds,]),
                       tot_umi = Matrix::colSums(seurat.object),
                       prop_mouse = n_mouse_umi / tot_umi,
                       prop_human = n_human_umi / tot_umi)

#Using our table we now determine the species of the cell in a new added column called "species"
#With case_when() we decide that if prop_mouse bigger that >95% the value in species is to be "mouse"
#if prop_human > 95% the value in species is to be "human", and if it's somewhere inbetween (e.g. 80% human, 20% mouse) the value is "mixed"

cell_species <- cell_species %>% mutate(species = case_when(prop_mouse > 0.95 ~ "mouse",prop_human > 0.95 ~ "human",TRUE ~ "mixed"))

#We can show the distribution of umis in a scatter plot with the color being the determined species
ggplot(cell_species, aes(n_human_umi, n_mouse_umi, color = species))+
      geom_point()+labs(title = "Human vs. Mouse UMIs")+theme_bw()

#We then add this information as metadata to our object
seurat.object <- AddMetaData(seurat.object, metadata = cell_species$species, col.name = "species")








#HIV transcriptional positive cell identification####
#This step basically follows the same principle as the species identification step

#We define every gene that has the prefix "mm10" or "GRCh38" or are the lenti virus genes "dsRed" or "EGFP" or "WPRE" as "HIV-" and the remaining as "HIV+"
#You can also flip this around and name all your HIV genes that you have defined in your reference and mark them as HIV+ and all remaining genes as HIV-
#!!!Check if the genes you are selecting using str_detect are actually the ones you want to select, some gene names might be similar and then falsely selected!!
#use str_subset(rownames(seurat.object), pattern = c("GENES OF INTEREST"))
#Sidenote: Adding * to the end of the name indicates it as part of a pattern, adding ^ before the name (e.g. "^gag-pol") indicates that there are no characters in front

gene_HIV <- ifelse(str_detect(rownames(seurat.object),"mm10*|GRCh38*|dsRed|EGFP|WPRE"),"HIV-", "HIV+")
hivpos_inds <- gene_HIV == "HIV+"
hivneg_inds <- gene_HIV == "HIV-"

#Here we again create a table by calcualting the percentages of HIV+/HIV- genes per cell
cell_hivstat <- tibble(n_hivpos_umi = Matrix::colSums(seurat.object[hivpos_inds,]),
                         n_hivneg_umi = Matrix::colSums(seurat.object[hivneg_inds,]),
                         tot_umi = Matrix::colSums(seurat.object),
                         prop_hivpos = n_hivpos_umi / tot_umi,
                         prop_hivneg = n_hivneg_umi / tot_umi)

#Here it becomes indivdual again and of course depends on your experimental design/dataset/subjective view
#In our datasets we saw a wide range of HIV+ umis (0 to >1000), so we decided on seperate ranges of HIV+ transkript positivity
#The nwely created column is called "status"
#Cells with >100 HIV+ umis are "HIV+ very high", 11-100 is "HIV+ high", 2-10 is "HIV+ low", and anything 1 and < is "HIV-"  
cell_hivstat <- cell_hivstat %>% mutate(status = case_when(n_hivpos_umi > 100 ~ "HIV+ very high",
                                                           n_hivpos_umi > 10 & n_hivpos_umi <= 100 ~ "HIV+ high",
                                                           n_hivpos_umi <= 1 ~ "HIV-",
                                                           n_hivpos_umi > 1 & n_hivpos_umi <= 10 ~ "HIV+ low"))

#Again we can show this in a scatter plot HIV+ umis vs. HIV- umis 
ggplot(cell_hivstat, aes(n_hivpos_umi, n_hivneg_umi, color = status))+ 
  geom_point()+scale_x_log10()+scale_y_log10()+
  labs("HIV+ vs. HIV- UMIs")+
  theme_bw()

#And then finally in the last step we add the "status" as meta data for each cell in our seurat object
seurat.object <- AddMetaData(seurat.object, metadata = cell_hivstat$status, col.name = "status")









#Lenti-RG transcriptionally positive cell identification####
#This follows the same principal as the two previous steps

gene_dsRed <- ifelse(str_detect(rownames(seurat.object),c("dsRed|WPRE|EGFP")),"LentiRG", "other")

dsRed_inds <- gene_dsRed == "LentiRG"
nondsRed_inds <- gene_dsRed == "other"

cell_dsRedstat <- tibble(n_dsRed_umi = Matrix::colSums(seurat.object[dsRed_inds,]),
                         n_other_umi = Matrix::colSums(seurat.object[nondsRed_inds,]),
                         tot_umi = Matrix::colSums(seurat.object),
                         prop_dsRed = n_dsRed_umi / tot_umi,
                         prop_other = n_other_umi / tot_umi)
  
#Any cell that has UMI >=1 of the Lenti-RG genes we determine as "Lenti-RG" - The column is named "dsRed expression" (technically kinda wrong)

cell_dsRedstat <- cell_dsRedstat %>% mutate(dsRed.expr = case_when(n_dsRed_umi >= 1 ~ "LentiRG",TRUE ~ "No expression"))

#The corresponding scatterplot
ggplot(cell_dsRedstat, aes(n_dsRed_umi, n_other_umi, color = dsRed.expr))+
  geom_point()+
  scale_x_log10()+scale_y_log10()+
  labs(title = "LentiRG UMIs")+
  theme_bw()

#We add this info as meta data to our seurat object
seurat.list[[i]] <- AddMetaData(seurat.object, metadata = cell_dsRedstat$dsRed.expr, col.name = "LentiRG.expression")


#QC - Part 1####
#We base the 1st QC round on nFeature_RNA count and mt% thresholds
#
#Here, the thresholds you decided on when looking at the previous violin plots come into play

seurat.object <- subset(seurat.object,subset= nFeature_RNA > 200 & nFeature_RNA < "UPPER THRESHOLD" & mitoHuCH38Percent < "mt-gene% cutoff")

#You can look at the violin plots post-subsetting again to check your thresholds
VlnPlot(seurat.object, features=c("nCount_RNA","mitoMM10Percent","mitoHuCH38Percent","nFeature_RNA"))+patchwork::plot_annotation(title = "After 1st QC") 

#As well as the scatterplots
FeatureScatter(seurat.object, feature1 = "nCount_RNA", feature2 = "mitoHuCH38Percent")+labs(title = "After 1st QC - human %mt")
FeatureScatter(seurat.object, feature1 = "nCount_RNA", feature2 = "mitoMM10Percent")+labs(title = "After 1st QC - mouse %mt")
FeatureScatter(seurat.object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+labs(title = "After 1st QC - nFeature_RNA")


#Identifying Doublets - DoubletFinder####

#This workflow follows the basic workflow for this package as described on their github
#(https://github.com/chris-mcginnis-ucsf/DoubletFinder)
#This is based on no ground-thruth (Because no multiplexing was performed)
#Therefore the assumed Doublet formation rate, depends on your 10X protocol and the number of expected cells
#e.g we used 5% since the cell -recovery rate for Namita's datasets was ~50.0000 and 10X gives a ~0.8% for 1000 cells
#So here you need to put in the appropriate expected multiplet rate for your dataset

#We first need to normalize and run a basic PCA on our dataset
#It's also a good idea to make a copy of the seurat.object before running this part

#Make a save copy
seurat.object.original <- seurat.object

#
seurat.object <- SCTransform(seurat.object)
seurat.object <- RunPCA(seurat.object)
seurat.object <- FindNeighbors(object = seurat.object, dims = 1:30)
seurat.object <- FindClusters(seurat.object) #For this we will use the default resolution settings
seurat.object <- RunUMAP(seurat.object, dims = 1:30)
  
#pK Identification (no ground-truth)
#Since there's no multiplexing we do not have any ground-truth (e.g. CMO-barcode-doublets)

sweep.res.seurat.list <- paramSweep_v3(seurat.object, PCs = 1:ncol(seurat.object@reductions$pca) , sct = TRUE)
sweep.stats.seurat.list <- summarizeSweep(sweep.res.seurat.list, GT = FALSE)
bcmvn.sweep.stats.seurat.list <- find.pK(sweep.stats.seurat.list)
  
#The default generated plot is generally really unhelpful so we generate our own
ggplot(bcmvn.sweep.stats.seurat.list, aes(pK, BCmetric, group = 1))+
       geom_point()+geom_line()+labs(title = "SeuratObject - pKs")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
#This line whill then automatically select optimal pK with the highest BCmetric
pK <- bcmvn.sweep.stats.seurat.list %>%
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]])) #Taking the top pK value

# Homotypic Doublet Proportion Estimate 
#Based on the number of cells in the dataset
nExp_poi <- round(0.05*nrow(seurat.object@meta.data))  ##Assuming 5% doublet formation rate - This rate is technically individual for each datasets depending on the expected recovery rate during the 10X protocol
  
annotations <- seurat.object@meta.data$seurat_clusters #based on the clustering results
homotypic.prop <- modelHomotypic(annotations)           
  
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
# Run DoubletFinder with varying classification stringencies
#Doublet/Singlet identity will be assigned to each cell
seurat.object <- doubletFinder_v3(seurat.object,
                                            PCs = 1:30,
                                            pN = 0.25,
                                            pK = pK,
                                            nExp = nExp_poi.adj, 
                                            reuse.pANN = FALSE,
                                            sct = TRUE)

#Showing Doublets/Singlet in the basic UMAP created
#The name for this column always contains the number of estimated doublets + variables used by DoubletFinder,
#which means it is never consistent. With str_subset we can select the column automatically without having to look at the objcet and 
#putting in the pull name of the new meta.data column

DimPlot(seurat.object,
        group.by = dplyr::first(str_subset(colnames(seurat.object@meta.data),pattern = "DF.classifications")),
        label=TRUE,repel=10,label.size=3)+
  NoLegend()+
  labs(title = paste0(dplyr::first(str_subset(colnames(seurat.object@meta.data),pattern = "DF.classifications"))))

##How many Doublets were found?
print(paste0("The SeuratObject had ",nExp_poi.adj," estimated Doublets"))
  
#Rename the column to get uniform meta.data columns
  
Indx<-str_which(colnames(seurat.object@meta.data),pattern="DF.classifications")
  
colnames(seurat.object@meta.data)[Indx] <- "Doublet.Singlet"

#We don't need the normalization and PCA/UMAP anymore so we take the meta.data and add it to the copy of seurat.object
#we made before running DoubletFinder

Doublet.Singlet <- FetchData(seurat.object, vars = c('Doublet.Singlet'), slot="counts") #In this table the cell.barcodes are the rownames, those are matching to the copy and are used to match the meta.data to the right cell

seurat.object.original<-AddMetaData(seurat.object.original,Doublet.Singlet) #The automatic column name will be the column name of the table

#We can now safely remove the SeuratObject we don't need by overwriting them

seurat.object<-seurat.object.original

rm(seurat.object.original)

##This is the point we're we would remove the "Doublet" cells
#In our analysis we kept th expected Doublets in

seurat.object <- subset(seurat.object,subset = Doublet.Singlet == "Singlet")

#QC - Part 2#####
#Now we're removing those cells with species "mouse" and "mixed"

seurat.object <- subset(seurat.object,subset = species == "human")

#We can create this table to check the number of remaining cells after our initial QC steps
table<-seurat.object@meta.data%>%dplyr::count(species)


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


