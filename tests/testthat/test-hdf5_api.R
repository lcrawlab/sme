test_that("end-to-end hdf5 api test", {
  # given
  hdf5_file <- tempfile()
  ds <- "gxg/0"
  mask <- 1:10
  
  # when
  create_hdf5_file(hdf5_file)
  write_hdf5_dataset(hdf5_file, ds, mask)
  result <- read_hdf5_dataset(hdf5_file, ds)


  # then
  expect_equal(result, mask)
})