# Writing data to file
filename = "testbed_file.csv"

# If file doesn't exist, create it...
if(!file.exists(filename)){
  file.create(filename)
} else { #...and if it exists, append data to it
  # Get user entries
  user_input = readline(prompt="Enter line: ")
  write(user_input,
        file=filename,
        append=TRUE)
}