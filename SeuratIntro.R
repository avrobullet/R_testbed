# Load the PBMC dataset
message('Loading HG19 gene data as a Seurat object...')
pbmc.data <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19/")
#pbmc.data <- Read10X_h5(data.dir = "filtered_gene_bc_matrices/hg19/")

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, 
                           project = "pbmc3k",
                           min.cells = 3, 
                           min.features = 200)
pbmc[['percent.mt']] <- PercentageFeatureSet(pbmc, '^MT-') 

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3)

# Return genes with RNA counts greater than 200 and less then 2500
message('Filter gene data to only retain genes with ...')
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalize data using log
message('Normalize data...')
pbmc <- NormalizeData(pbmc, 
                      normalization.method = "LogNormalize", 
                      scale.factor = 10000)

#' AD: This next step is about scaling the data to be used for PCA, which I don't 
#'     think NanoString has done at all prior to running PCAs.
message('Scale data using linear transform...')
pbmc <- ScaleData(pbmc, rownames(pbmc))
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Run PCA
message('Run Principal Component Analysis (PCA)...')
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Run UMAP
message('Run Uniform Manifold Approximation and Projection (UMAP)...')
pbmc <- RunUMAP(pbmc, dims = 1:10)

# Run kNN Leiden Clustering
message('Run Nearest Neighbour and Leiden clusterting...')
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5, algorithm = 4) # 4 = leiden

print(DimPlot(pbmc, reduction = "umap"))