## To run all tests enter devtools::test()

test_that("Loads valid image object", {
  x = LoadImage("~/Documents/test_tifs/HAB135_WT_Tom20-CAT-dapi_tile-1.czi_-_WT_Tom20-CAT-dapi_tile-1_02.tif")
  expect_type(x, "double")
})
