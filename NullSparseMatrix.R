
convertMatrix <- function(mat){
  rownames(blank_matrix) <- paste0("gene", 1:size)
  colnames(blank_matrix) <- paste0("cell_", 1:size)
  
  return(as(mat, "dgCMatrix"))
}

createBlankMatrix <- function(x){
  size <- 100
  blank_matrix <- Matrix(c(1:size),
                         byrow = TRUE, nrow = size, ncol = size, sparse = TRUE)
  rownames(blank_matrix) <- paste0("gene", 1:size)
  colnames(blank_matrix) <- paste0("cell_", 1:size)
  
  return(as(blank_matrix, "dgCMatrix"))
}

createMatrix <- function(rows = 2,
                         cols = 2,
                         matrix_type = "dgCMatrix") {
  new_matrix <- Matrix(1:(rows*cols), byrow = TRUE, nrow = rows, ncol = cols, sparse = TRUE)
  rownames(new_matrix) <- paste0("gene", 1:rows)
  colnames(new_matrix) <- paste0("cell_", 1:cols)
  
  return(as(new_matrix, matrix_type))
}

createSeuratBlank <- function(x){
  sem <- CreateSeuratObject(createBlankMatrix())
  # Return specific blank
  if(x=="counts"){
    return(sem@assays$RNA@counts)
  } else if (x=="data") {
    return(sem@assays$RNA@data)
  }
  else if (x=="i") {
    return(sem@assays$RNA@data@i)
  }
}