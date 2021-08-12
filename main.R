library(Seurat)
library(dplyr)
library(patchwork)
library(matrixStats)

# Object Referencing/Pointing (hacking)
#source(paste0(getwd(),"/ObjectReferencing.R"))
my_stats <- function(mat,iterations){
  mat_sum <- colSums(mat)
  mat_avg <- mat_sum/iterations
  mean_diff <- 0
  print(mat)
  print(mat_avg)
  for(j in 1:iterations){
    # Sum and squared
    mean_diff <- mean_diff + (mat[j]-mat_avg)^2
  }
  # Divide by total per col size and square root
  mean_diff <- sqrt(mean_diff/iterations)
  print(mean_diff)
}

#Seurat Demo Pipeline
iterations <- 10
seuratdemo_metrics <- matrix(nrow=iterations, ncol=1)
for(i in 1:iterations) {
  # Run metrics and pipeline
  start_time <- Sys.time()
  source(paste0(getwd(),"/SeuratIntro.R"))
  end_time <- Sys.time()
  # Save metrics
  seuratdemo_metrics[i,] <- end_time - start_time
}
print(seuratdemo_metrics)
my_stats(seuratdemo_metrics, iterations)

