library(Seurat)
library(dplyr)
library(patchwork)
library(matrixStats)
library(aws.s3)
library(sjmisc)

# Object Referencing/Pointing (hacking)
#source(paste0(getwd(),"/ObjectReferencing.R"))
#source(paste0(getwd(),"/SeuratIntro.R"))
#source(paste0(getwd(),"/AWSAccess.R"))
# demo <- c()
# demo_runs <- c('Run1',
#                'Run2',
#                'Run2',
#                'Run3')
# checkDuplicates <- function(run){
#   if(length(demo) == 0){
#     demo <<- c(demo, run)
#   }
#   else {
#     duplicate_found = FALSE
#     for(i in seq_len(length(demo))){
#       if(run == demo[i]){
#         duplicate_found = TRUE
#         break
#       }
#     }
#     if(!duplicate_found){demo <<- c(demo,run)}
#   }
# }
# sapply(demo_runs, checkDuplicates)
# print(demo)
words <- as.list(strsplit('My_name_is_BOB','_')[[1]])
print(words[[length(words)]])
