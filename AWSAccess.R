Sys.setenv(
  "AWS_ACCESS_KEY_ID" = "AKIAZBJQHHET6VIL65FJ",
  "AWS_SECRET_ACCESS_KEY" = "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",
  "AWS_DEFAULT_REGION" = "us-west-2"
)

# Variables
config_main <- list(
  resultspath = "100-plex/NSCLC_HIO/"
)

# Get S3 object from NanoString bucket
nanostring_bucket <- aws.s3::get_bucket("nanostring-external-smi-commercial")
# Convert NanoString S3 object into a data.table to better make use of its keys
nanostring_bucket_table <- data.table::rbindlist(nanostring_bucket)
# Only keep keys that are needed for the experiment
for(key in nanostring_bucket_table[, Key]){
  if(str_contains(key, config_main$resultspath)){ #str_contains(key, "100-plex/NSCLC_HIO/")
    print(key)
  }
}
