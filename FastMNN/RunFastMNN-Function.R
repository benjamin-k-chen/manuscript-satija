#Modified RunFastMNN function

#Before we can run FastMNN we need to modify the function however, since this version does not
#function with SCT-normalized datasets

##https://github.com/satijalab/seurat/issues/5329 (The error that RunFastMNN will throw)
#this is how they recommend to fix it
#trace('RunFastMNN', edit = T)
#replace line 25-28
# objects.sce <- lapply(X = object.list, FUN = function(x,f) { if (DefaultAssay(x) == 'SCT') { x = subset(x = x, features = f) indx <- match(rownames(x@assays$SCT@counts),rownames(x@assays$SCT@scale.data)) x@assays$SCT@scale.data <- x@assays$SCT@scale.data[indx,] }else{ x = subset(x = x, features = f) } return(as.SingleCellExperiment(x)) }, f = features)

RunFastMNN<-function (object.list, assay = NULL, features = 2000, reduction.name = "mnn", 
                       reduction.key = "mnn_", reconstructed.assay = "mnn.reconstructed", 
                       verbose = TRUE, ...) 
{
  if (!all(sapply(X = object.list, FUN = inherits, what = "Seurat"))) {
    stop("'object.list' must be a list of Seurat objects", 
         call. = FALSE)
  }
  if (length(x = object.list) < 2) {
    stop("'object.list' must contain multiple Seurat objects for integration", 
         call. = FALSE)
  }
  assay <- assay %||% DefaultAssay(object = object.list[[1]])
  for (i in 1:length(x = object.list)) {
    DefaultAssay(object = object.list[[i]]) <- assay
  }
  if (is.numeric(x = features)) {
    if (verbose) {
      message(paste("Computing", features, "integration features"))
    }
    features <- SelectIntegrationFeatures(object.list = object.list, 
                                          nfeatures = features, assay = rep(assay, length(x = object.list)))
  }
  objects.sce <- lapply(X = object.list, FUN = function(x,f) {
    if (DefaultAssay(x) == 'SCT') { x = subset(x = x, features = f) 
    indx <- match(rownames(x@assays$SCT@counts),rownames(x@assays$SCT@scale.data)) 
    x@assays$SCT@scale.data <- x@assays$SCT@scale.data[indx,] }
    else{ x = subset(x = x, features = f) }
    return(as.SingleCellExperiment(x)) }, f = features)
  integrated <- merge(x = object.list[[1]], y = object.list[2:length(x = object.list)])
  out <- do.call(what = batchelor::fastMNN, args = c(objects.sce, 
                                                     list(...)))
  rownames(x = SingleCellExperiment::reducedDim(x = out)) <- colnames(x = integrated)
  colnames(x = SingleCellExperiment::reducedDim(x = out)) <- paste0(reduction.key, 
                                                                    1:ncol(x = SingleCellExperiment::reducedDim(x = out)))
  integrated[[reduction.name]] <- CreateDimReducObject(embeddings = SingleCellExperiment::reducedDim(x = out), 
                                                       loadings = as.matrix(SingleCellExperiment::rowData(x = out)), 
                                                       assay = DefaultAssay(object = integrated), key = reduction.key)
  integrated[[reconstructed.assay]] <- CreateAssayObject(data = as(object = SummarizedExperiment::assay(x = out), 
                                                                   Class = "sparseMatrix"), )
  VariableFeatures(object = integrated[[reconstructed.assay]]) <- features
  Tool(object = integrated) <- S4Vectors::metadata(x = out)
  integrated <- LogSeuratCommand(object = integrated)
  return(integrated)
}
