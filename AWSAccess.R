Sys.setenv(
  "AWS_ACCESS_KEY_ID" = "AKIAZBJQHHET6VIL65FJ",
  "AWS_SECRET_ACCESS_KEY" = "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",
  "AWS_DEFAULT_REGION" = "us-west-2"
)

# Variables
config_main <- list(
  resultspath = '100-plex/NSCLC_HIO/'
)
config_loading <- list(
  annotfile_location = 'FOV annotation.csv', #16CPA annotation.xlsx
  annotfile = NULL
)
# runs_found <- list(
#   folder_path <- c(), #config_main$resultspath
#   slidefolders <- c(), #run_parameters[[5]]
#   votedfolders <- c(), #run_parameters[[6]]
#   version <-c(), #version_runname <- as.list(strsplit(run_parameters[[6]],'_')[[1]])
#                  #version_runname[[length(version_runname)]]
#   Run_name <- c(), #version_runname[[1]]
#   Slide_name <- c(), #run_parameters[[5]]
#   ISH_concentration <- c(),
#   Dash <- c(),
#   tissue <- c(), # tissue <- as.list(strsplit(run_parameters[[3]],'_')[[1]])
#                  # tissue <- as.list(strsplit(tissue[[2]],'-')[[1]])
#                  # tissue[[1]]
#   Panel <- c()
# )
runs_found <- c()

nanostring_bucket <- aws.s3::get_bucket("nanostring-external-smi-commercial")
nanostring_bucket_table <- data.table::rbindlist(nanostring_bucket)
# For all keys...
for(key in nanostring_bucket_table[,Key]){
  # ...only focus on keys matching the user's designated results path...
  if(str_contains(key, config_main$resultspath)){
    #print(key)
    # ...to load the necessary FOV annotations
    if(str_contains(key, config_loading$annotfile_location)){
      message(paste0('Annotation File Found: ',key))
      config_loading$annotfile <- aws.s3::s3read_using(FUN = read.csv, 
                                                       bucket = nanostring_bucket, 
                                                       object = key)
    }
    # ...to load the necessary data files
    if(str_contains(key, 'Run')){
      run_parameters <- as.list(strsplit(key,'/')[[1]])
      # Identify all unique runs
      checkDuplicates <- function(run){
        if(length(runs_found) == 0){
          runs_found <<- c(runs_found, run)
        }
        else {
          duplicate_found = FALSE
          for(i in seq_len(length(runs_found))){
            if(run == runs_found[i]){
              duplicate_found = TRUE
              break
            }
          }
          if(!duplicate_found){runs_found <<- c(runs_found,run)}
        }
      }
      sapply(run_parameters[6], checkDuplicates)
    }
    #'gem <- runloadingmodule(config_loading,
    #'                        loading_results_filepath = "intermediate_results/loading_results.RData",
    #'                        loading_input_hash_filepath = "intermediate_results/hashes/loading_hash.RData"))
  }
}
print(runs_found)

