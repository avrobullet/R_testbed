# library(Seurat)
# library(SeuratDisk)
# library(SeuratObject)
# library(dplyr)
# library(patchwork)
# library(parallel)
# library(Giotto)
# library(Matrix)
# library(reader)
# library(rhdf5)
# library(stringr)
library(DelayedArray)
library(HDF5Array)
library(ExperimentHub)

# Variables
filename <- "Seurat_demo.h5"

#setwd("C:/Users/alexa/R_testbed")
# Object Referencing/Pointing (hacking)
#source(paste0(getwd(),"/ObjectReferencing.R"))
#source(paste0(getwd(),"/SeuratIntro.R"))
#source(paste0(getwd(),"/AWSAccess.R"))
#source(paste0(getwd(),"/WriteToFile.R"))
#source(paste0(getwd(),"/NullSparseMatrix.R"))
#source(paste0(getwd(),"/VerifyFileContent.R"))
#source(paste0(getwd(),"/AppendMatrices.R"))
#source(paste0(getwd(),"/StreamHDF5Data.R"))

hub <- ExperimentHub()
fname <- hub[["EH1040"]]
tenx <- HDF5Array(filepath = fname, name = "counts")
tenx_subset <- tenx[1:1000,1:100000]
# for(i in 1:1) {
#   tenx_subset <- cbind(tenx_subset,tenx[1:1000,])
# }
print(dim(tenx_subset))
tenx_subset <- tenx_subset +1L
genes <- paste0("gene", 1:dim(tenx_subset)[1])
cells <- paste0("cell", 1:dim(tenx_subset)[2])
rownames(tenx_subset) <- genes
colnames(tenx_subset) <- cells

# Start time 
start_time = proc.time()[[3]]

# Create a Seurat object
sem <- CreateSeuratObject(count = as(tenx_subset, "dgCMatrix"))
sem <- NormalizeData(sem)
sem <- FindVariableFeatures(sem, selection.method = "vst", nfeatures = length(rownames(sem)))
sem <- ScaleData(sem, rownames(sem), center = FALSE)
sem <- RunPCA(sem, features = VariableFeatures(object = sem))
sem <- FindNeighbors(sem, dims = 1:10)
sem <- FindClusters(sem, resolution = 0.5)
sem <- RunUMAP(sem, dims = 1:10)
print(DimPlot(sem, reduction = "umap"))

# End time
end_time = proc.time()[[3]]

### Metrics ------------------------------------------------------------------
message(paste0("Total time (secs): ", end_time - start_time))


