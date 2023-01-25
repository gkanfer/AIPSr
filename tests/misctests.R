img = LoadImage("~/Documents/test_tifs/HAB135_WT_Tom20-CAT-dapi_1.tif")
nseg = segmentNucleus(img, index=3, minmaxnorm=TRUE, int=2, filter_size=50, offset=0.001, opensize=3,
                      small_obj=30, use_watershed=FALSE, rm_outliers=FALSE, out_p=0.95, displaymasks=TRUE)
cseg = segmentCyto(img, nseg$seg, index=2, int=10, filter_size=100,
                   offset=0.01, size_smooth=19, opensize=7, largeobj=100000, minmaxnorm=TRUE)

cyto_norm=norm_ch(image=img, index=2, minmaxnorm=TRUE)
EBImage::display(cyto_norm)

EBImage::display(cseg$norm)
EBImage::display(cseg$seg)
EBImage::display(nseg$seg)

str(cseg$norm)

table_shape = EBImage::computeFeatures.shape(cseg$seg,cseg$norm)
table_moment = EBImage::computeFeatures.moment(cseg$seg,cseg$norm)
table_basic = EBImage::computeFeatures.basic(cseg$seg,cseg$norm)
table_test <- as.data.frame(cbind(table_basic,table_moment,table_shape))
table_test$predict <- "Positive"
typeof(table_test)
