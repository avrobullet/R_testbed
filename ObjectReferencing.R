### Check object referencing ---------------------------------------------------

# Create object pointer function
newPointer=function(inputValue){ 
  object=new.env(parent=globalenv()) 
  object$value=inputValue 
  class(object)='pointer'
  
  return(object) 
} 

# Change contents of a matrix
alterMatrix <- function(m, k){
  for(i in 1:100) {
    for(j in 1:100){
      m[i,j] <- k
    }
  }
  return(m)
}

# Create variables
mat_a = matrix(nrow=100,ncol=100)
mat_b = newPointer(mat_a)
all_matrices = list(mat_a,
                    mat_b)

# Check matrices prior to being processed
# print(mat_a)
# print(mat_b$value)

# Process matrices
mat_a <- alterMatrix(mat_a,1)
mat_b$value <- alterMatrix(mat_b$value,2)

# Check matrices prior to being processed
print(mat_a)
print(mat_b$value)

#'VERIFIED: mat_b needs to be altered as its own entity (it lives in another
#'          environment). In addition, all_matrices seems to no update itself
#'          with mat_a after it was processed. Using pointers would make the 
#'          world of difference to my R coding!

# Check object size
print(object.size(mat_a))
print(object.size(mat_b))
print(object.size(all_matrices))