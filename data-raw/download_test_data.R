#!/usr/env Rscript

# Downloads a test data file from OpenScPCA and places it in the test directory

setwd(here::here())
test_file_path <- testthat::test_path("data", "scpca_sce.rds")
url <- "https://openscpca-test-data-release-public-access.s3.us-east-2.amazonaws.com/test/SCPCP000001/SCPCS000001/SCPCL000001_processed.rds" # nolint
download.file(url, test_file_path, mode = "wb")
