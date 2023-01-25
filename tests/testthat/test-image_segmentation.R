## To run all tests enter devtools::test()

test_that("Nuclei segmentation returns >3 mask objects, no outlier detection", {
  img = LoadImage("~/Documents/test_tifs/HAB135_WT_Tom20-CAT-dapi_1.tif")
  nseg = segmentNucleus(img, index=3, minmaxnorm=TRUE,
                        int=1, filter_size=50, offset=0.001, opensize=3,
                        small_obj=30, use_watershed=FALSE, 
                        rm_outliers=FALSE, out_p=0.95, displaymasks=TRUE)
  nf = EBImage::computeFeatures.moment(nseg$seg)
  expect(length(nf)>3)
})

test_that("Nuclei segmentation-remove outliers", {
  img = LoadImage("~/Documents/test_tifs/HAB135_WT_Tom20-CAT-dapi_tile-1.czi_-_WT_Tom20-CAT-dapi_tile-1_02.tif")
  nseg = segmentNucleus(img, index=3, minmaxnorm=TRUE,
                        int=1, filter_size=50, offset=0.001, opensize=3,
                        small_obj=30, use_watershed=FALSE, 
                        rm_outliers=TRUE, out_p=0.95, displaymasks=TRUE)
  nf = EBImage::computeFeatures.moment(nseg$seg)
  expect(length(nf)>3)
})

test_that("Cell segmentation returns >3 mask objects", {
  img = LoadImage("~/Documents/test_tifs/HAB135_WT_Tom20-CAT-dapi_tile-1.czi_-_WT_Tom20-CAT-dapi_tile-1_02.tif")
  nseg = segmentNucleus(img, index=1, minmaxnorm=TRUE,
                        int=2, filter_size=25, offset=0.01, opensize=3,
                        small_obj=30, use_watershed=FALSE, distmap_value=2,
                        rm_outliers=FALSE, out_p=0.95, displaymasks=TRUE)
  cseg = segmentCyto(img, nseg$seg, index=3, int=10, filter_size=100,
                     offset=0.01, size_smooth=19, opensize=7, largeobj=100000, minmaxnorm=TRUE)
  cf = EBImage::computeFeatures.moment(cseg$seg)
  expect(length(cf)>3)
})
