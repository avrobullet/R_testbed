Sys.setenv(
  "AWS_ACCESS_KEY_ID" = "AKIAZBJQHHET6VIL65FJ",
  "AWS_SECRET_ACCESS_KEY" = "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",
  "AWS_DEFAULT_REGION" = "us-west-2"
)

nanostring_bucket <- aws.s3::get_bucket("nanostring-external-smi-commercial")
for(bucket_key in nanostring_bucket$Contents){
  if(bucket_key == '100-plex/'){
    print(nanostring_bucket$Contents[bucket_key])
  }
}
# some_data <- s3read_using(FUN = read.csv,
#                           object = '100-plex/NSCLC_HIO/Run90089_Lung09c-HIO-1M-D03-S1/20210423_160913/20210423_160913_S1/Run90089_Iter1_PCIter0_TSR0.5_DF1_v5/FOV11/FOV011_ProcessingTimeLog.csv',
#                           bucket = nanostring_bucket)