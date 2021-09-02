library(DelayedArray)
library(HDF5Array)
library(ExperimentHub)
library(SingleCellExperiment)

# Create HDF5 file
createSeuratHDF5File <- function(filename){
  rhdf5::h5createFile(filename)
}

# Save Seurat Assays into individual HDF5 file
saveSeuratAssayHDF5File <- function(mat_type, mat_filename, mat_name) {
  message(paste0("\nSaving Seurat ",mat_name," as an HDF5 file..."))
  # Create file
  createSeuratHDF5File(mat_filename)
  # Save data
  rhdf5::h5createGroup(mat_filename,mat_name)
  rhdf5::h5write(mat_type@x,mat_filename,paste0(mat_name,"/x"))
  rhdf5::h5write(mat_type@i,mat_filename,paste0(mat_name,"/i"))
  rhdf5::h5write(mat_type@p,mat_filename,paste0(mat_name,"/p"))
  rhdf5::h5write(mat_type@Dim,mat_filename,paste0(mat_name,"/Dim"))
  rhdf5::h5write(mat_type@Dimnames,mat_filename,paste0(mat_name,"/Dimnames"))
  rhdf5::h5write(mat_type@factors,mat_filename,paste0(mat_name,"/factors"))
}

# Load Seurat Assays back from individual HDF5 file
loadSeuratAssayHDF5File <- function(mat_filename, mat_name) {
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
  rhdf5::h5createGroup(mat_filename,mat_name)
  rhdf5::h5write(mat_type@assay.used,mat_filename,paste0(mat_name,"/assay.used"))
  rhdf5::h5write(mat_type@i,mat_filename,paste0(mat_name,"/i"))
  rhdf5::h5write(mat_type@p,mat_filename,paste0(mat_name,"/p"))
  rhdf5::h5write(mat_type@Dim,mat_filename,paste0(mat_name,"/Dim"))
  rhdf5::h5write(mat_type@Dimnames,mat_filename,paste0(mat_name,"/Dimnames"))
  rhdf5::h5write(mat_type@factors,mat_filename,paste0(mat_name,"/factors"))
  # Return Seurat Blank
  return(createSeuratBlank(mat_name))
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


