library(rhdf5)
library(DelayedArray)
library(HDF5Array)
library(profvis)
library(ExperimentHub)
library(htmlwidgets)
library(Seurat)
library(Giotto)
library(ggplot2)
library(cowplot)

### SOURCE FILES ---------------------------------------------------------------
source("ProCogia_QCandNormalization.R")
source("ProCogia_SeuratObjectManagement.R")

### VARIABLES ------------------------------------------------------------------
N_SAMPLES = 1L # controls sample sizing
n_principal_components = 20
norm_data <- simulateData(N_SAMPLES)
norm_counts <- normalize_as_delayedarray(norm_data$data)

### PIPELINE -------------------------------------------------------------------
# Start time 
start_time = proc.time()[[3]]

### 1: CREATE SEURAT OBJECT ----------------------------------------------------
#' Needed to create Seurat blank without it calling integrity errors. Also, 
#' the addition of the normalized data is added to the Seurat as well.
createSpecializedSeuratObject <- function(rerun = FALSE){
  if (!isTRUE(rerun)){
    message("Create Seurat Object...")
    cbmc <- CreateSeuratObject(counts = norm_data$data)
    cbmc@assays$RNA@counts <- saveSeuratHDF5File(cbmc@assays$RNA@counts,"Seurat_Assay_RNA_Counts.h5","counts")
    cbmc@assays$RNA@data <- saveSeuratHDF5File(cbmc@assays$RNA@data,"Seurat_Assay_RNA_Data.h5","data")
    cbmc <- SetAssayData(cbmc, assay = "RNA", slot = "scale.data", new.data = as.matrix(norm_counts)) 
    config_hdf5$counts_location <- paste0(config_hdf5$hdf5_location,"/Seurat_Assay_RNA_Counts.h5")
    config_hdf5$data_location <- paste0(config_hdf5$hdf5_location,"/Seurat_Assay_RNA_Data.h5")
  } else {
    message("\nCreate Seurat Object from existing HDF5 files...")
    cbmc@assays$RNA@counts <- loadSeuratHDF5File(config_hdf5$counts_location,"counts")
    cbmc@assays$RNA@data <- loadSeuratHDF5File(config_hdf5$data_location,"data")
  }
  # Return specialized Seurat object
  return(cbmc)
}

### 2: SCALE NORMALIZED DATA ---------------------------------------------------
# Replace with SingleCellExperiment scale data function with DelayedArray
# cbmc <- ScaleData(cbmc, do.center = FALSE)

### 3: RUN PCA ON SEURAT OBJECT ------------------------------------------------
# Remove scale.data once PCA is done
runSpecializedPCA <- function(rerun = FALSE, cbmc){
  if (!isTRUE(rerun)){
    cbmc <- RunPCA(cbmc, features =  norm_data$genes, npcs = n_principal_components)
    cbmc@assays$RNA@scale.data <- as(c(1:4), "matrix")
  } else {
    message("\nConduct PCA from saved Seurat-HDF5 data...")
    # Fix...
  }
  # Return specialized Seurat object
  return(cbmc)
}

### 4: RUN NEAREST NEIGHBOUR ON SEURAT OBJECT ----------------------------------
runSpecializedNearestNeighbour <- function(rerun = FALSE, cbmc){
  if (!isTRUE(rerun)){
    cbmc <- FindNeighbors(cbmc, dims = 1:n_principal_components, annoy.metric = 'cosine' )
    cbmc@graphs$RNA_nn <- saveSeuratHDF5File(cbmc@graphs$RNA_nn,"Seurat_Graphs_RNA_nn.h5","RNA_nn")
  } else {
    message("\nConduct Nearest Neighbour from saved Seurat-HDF5 data...")
    # Fix...
  }
  # Return specialized Seurat object
  return(cbmc)
}

### 5: RUN CLUSTERING ALGORITHM ON SEURAT OBJECT -------------------------------
runSpecializedClusteringAlgorithm <- function(rerun = FALSE, cbmc){
  if (!isTRUE(rerun)){
    cbmc <- FindClusters(cbmc, resolution = 0.2, verbose = FALSE, algorithm = 1)
    cbmc@graphs$RNA_snn <- saveSeuratHDF5File(cbmc@graphs$RNA_snn,"Seurat_Graphs_RNA_snn.h5","RNA_snn")
  } else {
    message("\nConduct Clustering Algorithm from saved Seurat-HDF5 data...")
    # Fix...
  }
  # Return specialized Seurat object
  return(cbmc)
}

### 6: RUN UMAP ON SEURAT OBJECT -----------------------------------------------
runSpecializedUMAP <- function(rerun = FALSE, cbmc){
  if (!isTRUE(rerun)){
    cbmc <- RunUMAP(cbmc, dims = 1:n_principal_components, umap.method = "umap-learn")
    # cbmc@reductions$umap <- saveSeuratHDF5File()
  } else {
    
  }
  # Return specialized Seurat object
  return(cbmc)
}


sem <-  createSpecializedSeuratObject(FALSE)
sem <-  runSpecializedPCA(FALSE, sem)
sem <-  runSpecializedNearestNeighbour(FALSE, sem)
sem <-  runSpecializedClusteringAlgorithm(FALSE, sem)
sem <-  runSpecializedUMAP(FALSE, sem)

### 7: SAVE IMAGES FROM SEURAT OBJECT ------------------------------------------
save_plot("results.png", DimPlot(sem, label = TRUE, raster = FALSE))
# End time in seconds
end_time = proc.time()[[3]]

### METRICS --------------------------------------------------------------------
message(paste0("Total time: ", end_time - start_time, " seconds"))
message(paste0("Seurat object size: ", format(object.size(sem), units = "auto")))