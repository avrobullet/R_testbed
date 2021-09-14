# Keeps track of data being saved as an HDF5 file
config_hdf5 <- list(
  hdf5_location = paste0(getwd(),"/Seurat_HDF5_Files")
)

# Create HDF5 file
createSeuratHDF5File <- function(filename){
  hdf5_location <- paste0(config_hdf5$hdf5_location,"/",filename)
  rhdf5::h5createFile(hdf5_location)
  return(hdf5_location)
}

# Save Seurat RNA Assays into individual HDF5 file
saveSeuratHDF5File <- function(mat_type, mat_filename, mat_name) {
  message(paste0("\nSaving Seurat ",mat_name," as an HDF5 file..."))
  # Create file
  hdf5_mat_location <- createSeuratHDF5File(mat_filename)
  # Save data
  tryCatch({rhdf5::h5createGroup(hdf5_mat_location,mat_name)
    rhdf5::h5write(mat_type@i,hdf5_mat_location,paste0(mat_name,"/i"))
    rhdf5::h5write(mat_type@p,hdf5_mat_location,paste0(mat_name,"/p"))
    rhdf5::h5write(mat_type@Dim,hdf5_mat_location,paste0(mat_name,"/Dim"))
    rhdf5::h5write(mat_type@Dimnames,hdf5_mat_location,paste0(mat_name,"/Dimnames"))
    rhdf5::h5write(mat_type@x,hdf5_mat_location,paste0(mat_name,"/x"))
    rhdf5::h5write(mat_type@factors,hdf5_mat_location,paste0(mat_name,"/factors"))
  },
  error = function(e) { 
    print(paste0(mat_name," dataset exists.")) 
  })
  # Return empty Seurat
  mat_type@i <- as.integer(c(1:4))
  mat_type@p <- as.integer(c(1:4))
  mat_type@x <- as.integer(c(1:4))
  return(mat_type)
}

# Load Seurat Assays back from individual HDF5 file
loadSeuratHDF5File <- function(mat_filename, mat_name) {
  message(paste0("\nLoading Seurat ",mat_name," from specified HDF5 file..."))
  # Create empty matrix
  temp_sem <- createSeuratBlank(mat_name)
  str(temp_sem)
  # Load data
  temp_sem@x <- as.numeric(rhdf5::h5read(mat_filename, paste0(mat_name,"/x")))
  temp_sem@i <- as.integer(rhdf5::h5read(mat_filename, paste0(mat_name,"/i")))
  temp_sem@p <- as.integer(rhdf5::h5read(mat_filename, paste0(mat_name,"/p")))
  temp_sem@Dim <- as.integer(rhdf5::h5read(mat_filename, paste0(mat_name,"/Dim")))
  temp_sem@Dimnames <- rhdf5::h5read(mat_filename, paste0(mat_name,"/Dimnames"))
  temp_sem@factors <- rhdf5::h5read(mat_filename, paste0(mat_name,"/factors"))
  # Return data
  return(temp_sem)
}

# Save Seurat Graphs into individual HDF5 file
saveSeuratGraphHDF5File <- function(mat_type, mat_filename, mat_name) {
  message(paste0("\nSaving Seurat ",mat_name," as an HDF5 file..."))
  # Create file
  createSeuratHDF5File(mat_filename)
  # Save data
  tryCatch({rhdf5::h5createGroup(hdf5_mat_location,mat_name)
    rhdf5::h5write(mat_type@assay.used,mat_filename,paste0(mat_name,"/assay.used"))
    rhdf5::h5write(mat_type@i,hdf5_mat_location,paste0(mat_name,"/i"))
    rhdf5::h5write(mat_type@p,hdf5_mat_location,paste0(mat_name,"/p"))
    rhdf5::h5write(mat_type@Dim,hdf5_mat_location,paste0(mat_name,"/Dim"))
    rhdf5::h5write(mat_type@Dimnames,hdf5_mat_location,paste0(mat_name,"/Dimnames"))
    rhdf5::h5write(mat_type@x,hdf5_mat_location,paste0(mat_name,"/x"))
    rhdf5::h5write(mat_type@factors,hdf5_mat_location,paste0(mat_name,"/factors"))
  },
  error = function(e) { 
    print(paste0(mat_name," dataset exists.")) 
  })
  # Return Seurat 
  mat_type@i <- as.integer(c(1:4))
  mat_type@p <- as.integer(c(1:4))
  mat_type@x <- as.integer(c(1:4))
  return(mat_type)
}

# Load Seurat Graphs back from individual HDF5 file
loadSeuratGraphHDF5File <- function(mat_filename, mat_name) {
  message(paste0("\nLoading Seurat ",mat_name," from specified HDF5 file..."))
  # Create empty matrix
  temp_sem <- createSeuratBlank(mat_name)
  str(temp_sem)
  # Load data
  temp_sem@x <- rhdf5::h5read(mat_filename, paste0(mat_name,"/assay.used"))
  temp_sem@i <- as.integer(rhdf5::h5read(mat_filename, paste0(mat_name,"/i")))
  temp_sem@p <- as.integer(rhdf5::h5read(mat_filename, paste0(mat_name,"/p")))
  temp_sem@Dim <- as.integer(rhdf5::h5read(mat_filename, paste0(mat_name,"/Dim")))
  temp_sem@Dimnames <- rhdf5::h5read(mat_filename, paste0(mat_name,"/Dimnames"))
  temp_sem@factors <- rhdf5::h5read(mat_filename, paste0(mat_name,"/factors"))
  # Return data
  return(temp_sem)
}

# Save Seurat Graphs into individual HDF5 file
saveSeuratReductionHDF5File <- function(mat_type, mat_filename, mat_name) {
  message(paste0("\nSaving Seurat ",mat_name," as an HDF5 file..."))
  # Create file
  createSeuratHDF5File(mat_filename)
  # Save data
  rhdf5::h5createGroup(mat_filename,mat_name)
  
  # Return Seurat Blank
  return(createSeuratBlank(mat_name))
}

# Load Seurat Graphs back from individual HDF5 file
loadSeuratReductionHDF5File <- function(mat_filename, mat_name) {
  message(paste0("\nLoading Seurat ",mat_name," from specified HDF5 file..."))
  # Create empty matrix
  temp_sem <- createSeuratBlank(mat_name)
  str(temp_sem)
  # Load data
  
  # Return data
  return(temp_sem)
}
