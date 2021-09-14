### FUNCTIONS ------------------------------------------------------------------
#'Nehad: Revamped data normalization and quality control that doesn't involved
#'       Ptolemy and Giotto calls; just simply dcgMatrices!
#'AD: Currently, normalization is done via Seurat.
runQualityControl <- function(gem,
                              min_cell_counts = 4){
  #'Nehad: This is my simpleton way of seperating them please feel free to clean up
  gem_counts  <- gem@expression$rna$raw #AD: I fix!
  gem_genes <- gem@feat_ID$rna #AD: I fix!
  
  #'Nehad: identify negative genes columns
  #neg_genes_flag <- as.logical( sapply(gem_genes, function(y) grepl("NegPrb", y)))
  #'Nehad: identify FalseCode Genes
  #false_genes_flag <- as.logical( sapply(gem_genes, function(y) grepl("FalseCode", y)))
  #'Nehad: identify all other genes (anything but NegProb and FalseCode)
  #non_control_genes <- !(neg_genes_flag | false_genes_flag)
  #'Nehad: identify negative genes columns
  neg_counts <- gem_counts[as.logical( sapply(gem_genes, function(y) grepl("NegPrb", y))),]
  #'Nehad: identify FalseCode Genes
  false_counts <- gem_counts[as.logical( sapply(gem_genes, function(y) grepl("FalseCode", y))),]
  #'Nehad: identify all other genes that haven't been labelled with NegProb or FalseCode
  raw_counts <- gem_counts[!(as.logical( sapply(gem_genes, function(y) grepl("NegPrb", y))) | as.logical( sapply(gem_genes, function(y) grepl("FalseCode", y)))),]
  #'Nehad: get gene names/labels
  genes <- gem_genes[!(as.logical( sapply(gem_genes, function(y) grepl("NegPrb", y))) | as.logical( sapply(gem_genes, function(y) grepl("FalseCode", y))))]
  
  #'Nehad: remove low level genes ( can use highly variable genes to do a simliar thing)
  #'#'AD: ensured that the mean count of the dcgMatrix datastructure is
  #'      measuring the actually values (denoted as 'p')
  flagged_feats = Matrix::rowMeans(raw_counts) > Matrix::mean(neg_counts@p)
  genes <- genes[flagged_feats]
  raw_counts <- raw_counts[flagged_feats,]
  #'Nehad: remove cells with low counts
  flagged_cells = Matrix::colSums(raw_counts) > min_cell_counts
  raw_counts <- raw_counts[,flagged_cells]
  #'Nehad: normalize using totalcounts
  norm_counts <- sweep(raw_counts, 2, Matrix::colSums(raw_counts), "/") * mean(Matrix::colSums(raw_counts))
  return(list("raw_counts" = raw_counts,
              "norm_counts" = norm_counts,
              "feats" = genes,
              "flagged_cells" = flagged_cells))
}
#'Nehad: Normalize the DelayedArray data (done on-disk)
normalize_as_delayedarray <- function(counts){
  message("Normalizing data on DelayedArray...")
  darray_raw_counts <- DelayedArray(Matrix::t(counts))
  mean <- DelayedArray::mean(DelayedArray::rowSums(darray_raw_counts))
  col_normalized_P <- darray_raw_counts / DelayedArray::rowSums(darray_raw_counts)
  col_normalized_P <- Matrix::t(col_normalized_P * mean)
  return(col_normalized_P)
}
#'Nehad: Flag cells and genes (done on-RAM)
flag_cells_and_genes <- function(raw_counts, neg_counts, feats, min_cell_counts = 4){
  message("Flagging cells and genes...")
  #remove low level genes ( can use highly variable genes to do a simliar thing)
  flagged_feats = Matrix::rowMeans(raw_counts) > Matrix::mean(neg_counts)
  #remove cells with low counts 
  flagged_cells = Matrix::colSums(raw_counts) > min_cell_counts
  raw_counts <- raw_counts[,flagged_cells]
  results <-  list("feats" = flagged_feats, "flagged_cells" = flagged_cells)
  return(results)
}
#'Nehad: Add noise to the data (done on-RAM)
add_noise <- function(x, noise_rate = 0.05){
  message("Adding noise...")
  noise <- matrix( runif(length(x)), ncol= ncol(x) )
  noise[noise>noise_rate] <- 0 
  noise[noise>0] <- 1
  noisy_x <- x + noise
  return(noisy_x)
}
#'Nehad: Data filtered and corrected provided by Zhi to create +1 million cells
simulateData <- function(N_SAMPLES){
  load(paste0(getwd(),"/ProCogia_Pipeline/CPA data.RData"))
  raw_list  = c() 
  neg_list = c()
  hdf5_file_location <- paste0(getwd(),"/Seurat_HDF5_Files")
  
  #some variables to initialize
  neg_list = c()
  
  #here they combined all the genes in one matrix (experiment + control)
  # This is my simpleton way of seperating them please feel free to clean up :)
  raw_gene_counts  <- GemiObj@expression$rna$raw
  neg_counts <- GemiObj@expression$negprobes$raw
  genes <- rownames(raw_gene_counts)
  rm(GemiObj)
  gc()
  
  for (i in 1:N_SAMPLES){
    #generate random sampling integers
    rand_rows <- sample(1:ncol(raw_gene_counts), ncol(raw_gene_counts), replace = FALSE)
    #use rand_rows to randomly sample raw_genes and neg_counts
    new_raw <- raw_gene_counts[,rand_rows]
    new_neg <- neg_counts[, rand_rows]
    
    #add noise
    new_raw <- add_noise(new_raw)
    new_neg <- add_noise(new_neg)
    
    #append matrix column wise
    raw_list <- cbind(raw_list, new_raw)
    neg_list <- cbind(neg_list, new_neg)
  }
  
  #lets rename cell IDs to avoid duplicate names
  colnames(raw_list) <- paste0("cell_", 1:ncol(raw_list),"")
  
  #min_cell_counts = 5 100-plex, 20 for 1k-plex
  res <- flag_cells_and_genes(raw_list, neg_list, genes, min_cell_counts = 20)
  raw_list <- raw_list[res$feats,res$flagged_cells]
  neg_list <- neg_list[,res$flagged_cells] #'AD: this isn't used anywhere
  rm(new_neg)
  rm(new_raw)
  gc()
  return(list("data" = raw_list,
              "genes"=genes[res$feats]))
}
