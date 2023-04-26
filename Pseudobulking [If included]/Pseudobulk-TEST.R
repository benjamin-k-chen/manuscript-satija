#Pseudobulking - TEST Skript
#
#workflow adapted from: https://github.com/kpatel427/YouTubeTutorials/blob/main/singleCell_pseudoBulk.R
#
#We will be using the raw counts form the datasets after the final workflow
#
#PDF FILE
#
#FILEPATH SEURATLIST
#
#
library(Seurat)
library(SeuratObject)
library(dplyr)
library(MAST)#For Later 3DE expression analysis
library(scCustomize)
library(tidyverse)
library(Libra)
#
#Load seurat.list
seurat.list.merged<-merge(seurat.list[[3]],c(seurat.list[[7]],seurat.list[[8]]))

#What do we want to compare across the mice?
#GFP vs. dsRed
#HIV+ high vs. HIV-
#
#
#
#
# pseudo-bulk workflow
# Acquiring necessary metrics for aggregation across cells in a sample
# 1. counts matrix - sample level
# counts aggregate to sample level

#We could combine the meta.data colums MouseID and Condition in one column
#In the tutorial they group by cell tpye (-> predicted.celltype.l2) and sample (contains sampleID+Condition)
#We want to look at HIV+ across three acute mice?

cts <- AggregateExpression(seurat.list.merged, 
                           group.by = c("status","MouseID"),
                           assays = 'RNA',
                           slot = "counts",
                           return.seurat = FALSE)

cts <- cts$RNA

# transpose
cts.t <- t(cts)


# convert to data.frame
cts.t <- as.data.frame(cts.t)

#But we only want to compare the HIV+ high across the acute
#Split by HIV+ high

# get values where to split
splitRows <- gsub('_.*', '', rownames(cts.t))


# split data.frame
cts.split <- split.data.frame(cts.t,
                              f = factor(splitRows))


# fix colnames and transpose

cts.split.modified <- lapply(cts.split, function(x){
  rownames(x) <- gsub('.*_(.*)', '\\1', rownames(x))
  t(x)
  
})

# Let's run DE analysis with only HIV+ high cells of A1,A2,A3
# 1. Get counts matrix
counts_hivhigh <- cts.split.modified$`HIV+ high`

# 2. generate sample level metadata
colData <- data.frame(MouseID = colnames(counts_hivhigh))

colData <- colData %>%
  mutate(Condition = 'Acute 15dpi')%>%
  column_to_rownames(var = 'MouseID')

#3.Create a SeuratObject
testseurat<-CreateSeuratObject(counts_hivhigh,meta.data = colData)

markers<-FindAllMarkers(testseurat)
