#'AD: A piggy-back Seurat object that just holds the processed data.
#'@param mat: The matrix being created
#'@param row_amount: The matrix's user selected rows
#'@param col_amount: The matrix's user selected cols
dummySeurat <- function(row_amount,col_amount){
  mat <- matrix(c(1:row_amount*col_amount),nrow=row_amount,ncol=col_amount)
  
  rownames(mat) <- paste0("gene_", 1:dim(mat)[1])
  colnames(mat) <- paste0("cell_", 1:dim(mat)[1])
  return(as(mat, "dgCMatrix"))
}
dummySeurat(2,2)