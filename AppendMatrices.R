# Concatenate a matrix over and over again
mat1 <- matrix(1:12, ncol = 3)                    # Create first example matrix
mat2 <- matrix(21:35, ncol = 3)                   # Create second example matrix

# Row bind the matrices together (add mat2 as extra rows for mat1)
mat_combined1 <- rbind(mat1, mat2)                # Apply rbind function
print(mat_combined1)                                     # Print concatenated matrix