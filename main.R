library(Seurat)
library(dplyr)
library(patchwork)
library(hash)
library(matrixStats)

# Object Referencing/Pointing (hacking)
#source(paste0(getwd(),"/ObjectReferencing.R"))
# Seurat Demo Pipeline
# iterations <- 1
# seuratdemo_metrics <- matrix(nrow=iterations, ncol=1)
# for(i in 1:iterations) {
#   # Run metrics and pipeline
#   start_time <- Sys.time()
#   source(paste0(getwd(),"/SeuratIntro.R"))
#   end_time <- Sys.time()
#   # Save metrics
#   seuratdemo_metrics[i,] <- end_time - start_time
# }
# print(seuratdemo_metrics)
# avg_score <- sum(seuratdemo_metrics)
# # avg_score = colMean(seuratdemo_metrics)
# # print(sapply(seuratdemo_metrics, sd))
# hist(seuratdemo_metrics)
x <- matrix(rep(1:9), 3, 3)
print(colSums(x)/3)