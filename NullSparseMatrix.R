library(Seurat)

convertMatrix <- function(mat){
  rownames(blank_matrix) <- paste0("gene", 1:size)
  colnames(blank_matrix) <- paste0("cell_", 1:size)
  
  return(as(mat, "dgCMatrix"))
}

createBlankMatrix <- function(x){
  size <- 100
  blank_matrix <- Matrix(c(1:size),
                         byrow = TRUE, nrow = size, ncol = size, sparse = TRUE)
  rownames(blank_matrix) <- paste0("gene", 1:size)
  colnames(blank_matrix) <- paste0("cell_", 1:size)
  
  return(as(blank_matrix, "dgCMatrix"))
}

createMatrix <- function(rows = 2,
                         cols = 2,
                         matrix_type = "dgCMatrix") {
  new_matrix <- matrix(1:(rows*cols), byrow = TRUE, nrow = rows, ncol = cols)
  rownames(new_matrix) <- paste0("gene", 1:rows)
  colnames(new_matrix) <- paste0("cell_", 1:cols)
  
  return(as(new_matrix, matrix_type))
}

createSeuratBlank <- function(x){
  nrow <- 1000
  ncol <- 1000
  dummy_counts <- createMatrix(nrow,ncol)
  metadata <- list(rows = rownames(dummy_counts),
                   cols = colnames(dummy_counts))
  print(metadata$rows)
  sem <- CreateSeuratObject(createMatrix())
  #'AD: Somehow, the data slot is the comparator for adding datat to the RNA slot
  sem@assays$RNA@data@Dim <- as.integer(c(length(rownames(dummy_counts)),length(colnames(dummy_counts))))
  sem@assays$RNA@data@Dimnames <- list(rownames(dummy_counts),colnames(dummy_counts))
  sem <- SetAssayData(sem, assay = "RNA", slot = "scale.data", new.data = as.matrix(dummy_counts))
  
  # sem <- FindVariableFeatures(sem, selection.method = "vst", nfeatures = length(rownames(sem@assays$RNA@data@Dimnames[2])))
  # sem <- ScaleData(sem, assay = "scale.data", rownames(sem), center = FALSE)
  # sem <- RunPCA(sem, assay = "RNA",  slot = "scale.data")
  # sem <- FindNeighbors(sem, dims = 1:10)
  # sem <- FindClusters(sem, resolution = 0.5)
  # sem <- RunUMAP(sem, dims = 1:10)
  # print(DimPlot(sem, reduction = "umap"))
  
}
str(createMatrix())