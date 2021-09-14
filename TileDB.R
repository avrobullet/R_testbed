message("Running TileDB.R...")

# Create dimension
dim1 <- tiledb_dim("dim1", c(1L, 4L), 2L, "INT32")
# String dimenions: no values for domain and extent
dim2 <- tiledb_dim("dim2", NULL, NULL, "ASCII")


# Create domain with two dimensions
dom <- tiledb_domain(dims = c(dim1, dim2))