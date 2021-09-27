#'AD: Write performance stats to a dataframe

usedMemory <- function(x){
  used_memory <- total_memory - as.numeric(system("awk '/MemFree/ {print $2}' /proc/meminfo",intern=TRUE))
  return(used_memory)
}
currentTime <- function(x){
  return(proc.time()[[3]])
}

stats <- NULL
createHandler <- function(x){
  if(is.null(stats)){
    message("Create performance data handler...")
    stats <- data.frame(time=c(""), memory=c(""))
  } else {
    message("Collect performance data...")
  }
  return(stats)
}

