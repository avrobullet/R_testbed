# Sys.setenv(
#   "AWS_ACCESS_KEY_ID" = "AKIAZBJQHHET6VIL65FJ",
#   "AWS_SECRET_ACCESS_KEY" = "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",
#   "AWS_DEFAULT_REGION" = "us-west-2"
# )
# 
# # Variables
# config_main <- list(
#   resultspath = '100-plex/NSCLC_HIO/'
# )
# config_loading <- list(
#   annotfile_location = 'FOV annotation.csv', #16CPA annotation.xlsx?
#   annotfile = NULL
# )
# runs_found <- c()
# # Identify all unique runs
# checkDuplicates <- function(run){
#   if(length(runs_found) == 0){
#     runs_found <<- c(runs_found, run)
#   }
#   else {
#     duplicate_found = FALSE
#     for(i in seq_len(length(runs_found))){
#       if(run == runs_found[i]){
#         duplicate_found = TRUE
#         break
#       }
#     }
#     if(!duplicate_found){runs_found <<- c(runs_found,run)}
#   }
# }
# 
# nanostring_bucket <- aws.s3::get_bucket("nanostring-external-smi-commercial")
# nanostring_bucket_table <- data.table::rbindlist(get_bucket(bucket = nanostring_bucket))
# # tempfile <- tempfile()
# # save_object(object = 's3://nanostring-external-smi-commercial/100-plex/Run90070_100pIO-1nM-D3_20210317_155221_S1/CellStatsDir/Run90070_20210317_155221_Cell_Stats_F001.csv',
# #             file = tempfile)
# # write.csv(tempfile)
# # Convert to mlapply (or some other multithreading function)
# # For all keys...
# for(key in nanostring_bucket_table[,Key]){
#   # ...only focus on keys matching the user's designated results path...
#   if(str_contains(key, config_main$resultspath)){
#     # ...to load the necessary FOV annotations.
#     if(str_contains(key, config_loading$annotfile_location)){
#       message(paste0('Annotation File Found: ',key))
#       config_loading$annotfile <- aws.s3::s3read_using(FUN = read.csv,
#                                                        bucket = nanostring_bucket,
#                                                        object = key)
#     }
#     # ...to load the necessary data files.
#     if(str_contains(key, 'Run')){
#       # Break down Run# keys
#       run_parameters <- as.list(strsplit(key,'/')[[1]])
#       # print(key)
#       # Identify Run# experiment and create corresponding directories if they don't exist
#       # Create main experiment folder (Ex: Run####)
#       if(!file.exists(run_parameters[[3]])){
#         dir.create(run_parameters[[3]])
#       }
#       # Cell stats directory (CellStatsDir)
#       if(str_contains(key,"/CellStatsDir")){
#         if(!file.exists(paste0(run_parameters[[3]],'/CellStatsDir'))){
#           dir.create(paste0(run_parameters[[3]],'/CellStatsDir'))
#         } else { # Download data if directory exists
#           # Save old (previous) directory
#           oldwd <- getwd()
#           setwd(paste0(getwd(),'/',run_parameters[[3]],'/CellStatsDir'))
#           # If data is .csv
#           if(str_contains(run_parameters[[7]],".csv")){
#             save_object(key, file = run_parameters[[7]], bucket = nanostring_bucket)
#           }
#           # Reset to previous directory
#           setwd(oldwd)
#         }
#       }
#       # Cell Gene expression directory (Run####_Iter#_PCIter#_TSR_v#)
#       if(grepl('Run[0-9]+_Iter[0-9]+_PCIter[0]+_TSR[0-9]+.[0-9]+_DF[0-9]+_v[0-9]',run_parameters[6])){
#         if(!file.exists(paste0(run_parameters[[3]],'/',run_parameters[6]))){
#           dir.create(paste0(run_parameters[[3]],'/',run_parameters[6]))
#         } else if(!file.exists(paste0(run_parameters[[3]],'/',run_parameters[6],'/',run_parameters[7]))){
#           dir.create(paste0(run_parameters[[3]],'/',run_parameters[6],'/',run_parameters[7]))
#         } else if(!file.exists(paste0(run_parameters[[3]],'/',run_parameters[6],'/',run_parameters[7],'/',run_parameters[8]))){
#           dir.create(paste0(run_parameters[[3]],'/',run_parameters[6],'/',run_parameters[7],'/',run_parameters[8]))
#         } else { # Download data if directories exists
#           # Save old (previous) directory
#           oldwd <- getwd()
#           setwd(paste0(getwd(),'/',run_parameters[6],'/',run_parameters[7],'/',run_parameters[8]))
#           # If data is .csv
#           if(str_contains(run_parameters[[9]],".csv")){
#             save_object(key, file = run_parameters[[9]], bucket = nanostring_bucket)
#           }
#           # Reset to previous directory
#           setwd(oldwd)
#         }
#       }
#     }
#   }
# }
# 
