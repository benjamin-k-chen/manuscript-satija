#Function - RenameGenes
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
