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
source("ProCogia_Autoencoder.R")

### VARIABLES ------------------------------------------------------------------
N_SAMPLES = 4L #250L fo # controls sample sizing
total_memory <- as.numeric(system("awk '/MemTotal/ {print $2}' /proc/meminfo",intern=TRUE))
n_principal_components = 20
norm_data <- simulateData(N_SAMPLES)
norm_counts <- normalize_as_delayedarray(norm_data$data)
message(paste0("Simulated data size: ", ncol(norm_counts)))
start_time <- proc.time()[[3]]

### FUNCTIONS  ------------------------------------------------------------------
usedMemory <- function(x){
  used_memory <- total_memory - as.numeric(system("awk '/MemFree/ {print $2}' /proc/meminfo",intern=TRUE))
  return(used_memory)
}
currentTime <- function(x){
  return(proc.time()[[3]])
}

### 1: CREATE SEURAT OBJECT ----------------------------------------------------
#' Needed to create Seurat blank without it calling integrity errors. Also, 
#' the addition of the normalized data is added to the Seurat as well.
createSpecializedSeuratObject <- function(rerun = FALSE){
  if (!isTRUE(rerun)){
    message("Create Seurat Object...")
    cbmc <- CreateSeuratObject(counts = norm_data$data)
    cbmc@assays$RNA@counts <- saveSeuratHDF5File(cbmc@assays$RNA@counts,"counts","Seurat_Assay_RNA_Counts.h5")
    cbmc@assays$RNA@data <- saveSeuratHDF5File(cbmc@assays$RNA@data,"data","Seurat_Assay_RNA_Data.h5")
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
    cbmc <- approximatePCA(cbmc)
    cbmc@assays$RNA@scale.data <- as(c(1:4), "matrix")
    config_hdf5$PCA_location <- paste0(config_hdf5$hdf5_location,"/Seurat_Reduc_Autoencoder.h5")
  } else {
    message("\nConduct PCA from saved Seurat-HDF5 data...")
    cbmc@assays$RNA@counts <- loadSeuratHDF5File(config_hdf5$counts_location,"counts")
    # Fix...
  }
  # Return specialized Seurat object
  return(cbmc)
}

### 4: RUN NEAREST NEIGHBOUR ON SEURAT OBJECT ----------------------------------
runSpecializedNearestNeighbour <- function(rerun = FALSE, cbmc){
  if (!isTRUE(rerun)){
    cbmc <- FindNeighbors(cbmc, reduction = "autoencoder", dims = 1:n_principal_components, annoy.metric = 'cosine' )
    cbmc@graphs$RNA_nn <- saveSeuratHDF5File(cbmc@graphs$RNA_nn,"RNA_nn","Seurat_Graphs_RNA_nn.h5")
    config_hdf5$RNA_nn_location <- paste0(config_hdf5$hdf5_location,"/Seurat_Graphs_RNA_nn.h5")
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
    cbmc <- FindClusters(cbmc, reduction = "autoencoder", resolution = 0.2, verbose = FALSE, algorithm = 1)
    cbmc@graphs$RNA_snn <- saveSeuratHDF5File(cbmc@graphs$RNA_snn,"RNA_snn","Seurat_Graphs_RNA_snn.h5")
    config_hdf5$RNA_snn_location <- paste0(config_hdf5$hdf5_location,"/Seurat_Graphs_RNA_snn.h5")
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
    cbmc <- approximateUMAP(cbmc)
    #cbmc <- RunUMAP(cbmc, dims = 1:n_principal_components, umap.method = "umap-learn")
    # cbmc@reductions$umap <- saveSeuratHDF5File()
    config_hdf5$UMAP_location <- paste0(config_hdf5$hdf5_location,"/Seurat_Graphs_UMAP.h5")
  } else {
    
  }
  # Return specialized Seurat object
  return(cbmc)
}

sem <-  createSpecializedSeuratObject(FALSE)
message(paste0("Post-Seurat Object Time Point: ",currentTime()," secs"))
message(paste0("Post-Seurat Object RAM State: ",usedMemory()," kb"))
sem <-  runSpecializedPCA(FALSE, sem)
message(paste0("Post-PCA Object Time Point: ",currentTime()," secs"))
message(paste0("Post-PCA RAM State: ",usedMemory()," kb"))
sem <-  runSpecializedNearestNeighbour(FALSE, sem)
message(paste0("Post-NN Object Time Point: ",currentTime()," secs"))
message(paste0("Post-NN RAM State: ",usedMemory()," kb"))
sem <-  runSpecializedClusteringAlgorithm(FALSE, sem)
message(paste0("Post-Cluster Object Time Point: ",currentTime()," secs"))
message(paste0("Post-Cluster RAM State: ",usedMemory()," kb"))
sem <-  runSpecializedUMAP(FALSE, sem)
message(paste0("Post-UMAP Object Time Point: ",currentTime()," secs"))
message(paste0("Post-UMAP RAM State: ",usedMemory()," kb"))

### 7: SAVE IMAGES FROM SEURAT OBJECT ------------------------------------------
save_plot("results.png", DimPlot(sem, reduction = "approxumap", label = TRUE, raster = FALSE))
# End time in seconds
end_time <- proc.time()[[3]]

### METRICS --------------------------------------------------------------------
message(paste0("Total time: ", end_time - start_time, " seconds"))
message(paste0("Seurat object size: ", format(object.size(sem), units = "auto")))