test_that("Extract features works", {
  img = LoadImage("~/Documents/test_tifs/HAB135_WT_Tom20-CAT-dapi_1.tif")
  nseg = segmentNucleus(img, index=3, minmaxnorm=TRUE,
                        int=2, filter_size=25, offset=0.01, opensize=3,
                        small_obj=30, use_watershed=FALSE, distmap_value=2,
                        rm_outliers=FALSE, out_p=0.95, displaymasks=TRUE)
  cseg = segmentCyto(img, nseg$seg, index=2, int=10, filter_size=100,
                     offset=0.01, size_smooth=19, opensize=7, largeobj=100000, minmaxnorm=TRUE)
  table = extractFeatures(img, cseg$seg)
  expect(nrow(table$Ts.mix)>3)
})

# test_that("Cell selection/classification runs without error", {
#   img = LoadImage("~/Documents/test_tifs/HAB135_WT_Tom20-CAT-dapi_1.tif")
#   nseg = segmentNucleus(img, index=1, minmaxnorm=TRUE,
#                         int=2, filter_size=25, offset=0.01, opensize=3,
#                         small_obj=30, use_watershed=FALSE, distmap_value=2,
#                         rm_outliers=FALSE, out_p=0.95, displaymasks=TRUE)
#   cseg = segmentCyto(img, nseg$seg, index=3, int=10, filter_size=100,
#                      offset=0.01, size_smooth=19, opensize=7, largeobj=100000, minmaxnorm=TRUE)
#   table = extractFeatures(img, cseg$seg)
#   expect_output(pickCells(nseg$seg, cseg$norm, nseg$features, table$Ts.mix, int=10, font_size=0.7, label_class="Positive", display_select=TRUE))
# })

test_that("Pick cells correct output, selected cells", {
    img = LoadImage("~/Documents/test_tifs/HAB135_WT_Tom20-CAT-dapi_1.tif")
    nseg = segmentNucleus(img, index=3, minmaxnorm=TRUE,
                          int=2, filter_size=25, offset=0.001, opensize=3,
                          small_obj=30, use_watershed=FALSE, distmap_value=2,
                          rm_outliers=FALSE, out_p=0.95, displaymasks=TRUE)
    cseg = segmentCyto(img, nseg$seg, index=2, int=10, filter_size=100, offset=0.01, size_smooth=19, opensize=7, largeobj=100000, minmaxnorm=TRUE)
    table = extractFeatures(img, cseg$seg)
    out=pickCells(nseg$seg, cseg$norm, nseg$features, table$Ts.mix, int=10, font_size=0.7, label_class="Postive", display_select=TRUE)
    expect(nrow(out$table_train)>3)
})

#! need demo classification data
# test_that("Test SVM model function builds model & returns its accuracy from training data", {
#   
#   total_table = rbind(pos, neg)
# })