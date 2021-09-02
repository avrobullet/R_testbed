#'AD: Check FOV .csv file contents. If there aren't anymore row indices greater
#'    than one (1), simply skip the FOV folder from having it be loaded into the 
#'    pipeline.
#'@param filedir: directory of the file location
#'@param filename: name of every file in the file location
#'@param filename_path: path to the file itself to be read 
checkFileSize <- function(filedir){
  for(filename in list.files(path = filedir)) {
    filename_path <- paste0(getwd(),"/",filedir,"/",filename) 
    if(length(n.readLines(filename_path,n=5, header=TRUE)) <= 1){
      message(paste0("ATTENTION: ", filename," is not fully populated. Do not add associated data for analysis"))
      next
    }
  }
}