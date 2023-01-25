#' 1.1 Nucleus Segmentation
#' @param image EBImage object
#' @param index select which channel (1=ch1, 2=ch2..) integer value
#' @param minmaxnorm whether or not to normalize each channel by min/max intensity (TRUE or FALSE)
#' @param int intensity of image (integer)
#' @param filter_size size (x and y) of threshold window
#' @param offset ?
#' @param opensize ?
#' @param small_obj Integer, all objects smaller than value are removed 
#' @param use_watershed TRUE/FALSE, whether to use EBImage watershed for segmenting, instead of simpler bwlabel function
#' @param distmap_value If watershed=TRUE, edge sensitivity of segmentation
#' @param rm_outliers TRUE/FALSE ?
#' @param out_p ?
#' @return list: segmented nuclei mask, xy features, (display masks on image)
#' @export
segmentNucleus <- function(image, index=1, minmaxnorm=TRUE,
                           int=2, filter_size=16, offset=0.04, opensize=3,
                           small_obj=30, use_watershed=FALSE, distmap_value=2,
                           rm_outliers=TRUE, out_p=0.95, displaymasks=TRUE) {
  norm_nuc <- norm_ch(image=image, index=index, minmaxnorm=minmaxnorm)
  mask_nuc <- mask_ch(norm_nuc, filter_size, offset, opensize, shape="diamond")
  if (use_watershed) {
    seg_nuc = EBImage::watershed(distmap(mask_nuc), distmap_value)
  } else {
    seg_nuc = EBImage::bwlabel(mask_nuc)
  }
  chkpt(seg_nuc)
  #Remove small objects
  nf=EBImage::computeFeatures.shape(seg_nuc)
  nr=which(nf[,2]<small_obj)
  nseg=EBImage::rmObjects(seg_nuc, nr)
  chkpt(nseg)
  #Remove Outliers
  if (rm_outliers) {
    int.dapi<-EBImage::computeFeatures.basic(nseg, image)
    y=which(outliers::scores(int.dapi[,1], type="z", prob=0.95))
    tp=as.numeric(attr(y, "names"))
    if (length(tp)<1) {
      stop("Outlier detection failed")
    }
    nseg<-EBImage::rmObjects(nseg,tp)
  }
  dfshape<-EBImage::computeFeatures.shape(nseg)
  xy<-EBImage::computeFeatures.moment(nseg)[,c('m.cx', 'm.cy')]
  ##chkpt2(xy)

  if (displaymasks){
    seg_CH1mask2 = EBImage::paintObjects(mask_nuc,EBImage::toRGB(norm_nuc),opac=c(1, 1),col=c("red",NA),thick=TRUE,closed=TRUE)
    seg_mask<-EBImage::paintObjects(nseg,EBImage::toRGB(norm_nuc),opac=c(1, 1),col=c("red",NA),thick=TRUE,closed=TRUE)
  }
  return(list(seg=nseg, features=xy, mask_start=seg_CH1mask2, mask_final=seg_mask))
}

#' 1.2 Cell Segmentation
#' Segments cells using nuclei seeds and cytosolic marker
#' @param x loaded image
#' @param y nseg returned from segmentNucleus
#' @param index Integer. Which channel (1,2...) has cytosolic / phenotype image
#' @param largeobj Integer. Maximum size (area) of cell included in final segmentation.
#' @inheritParams norm_ch
#' @return List. Nucleus segmentation mask (nuc_seg), cytosolic segmentation mask (seg), mask features for nucleus segmentation (xy_nuc) and cytosolic segmentation (xy_cseg), normalized single cytosolic channel image (norm)
#' @export

segmentCyto <- function(x, y, index=2, int=40, filter_size=10, offset=0.1, size_smooth=19,
                   opensize=7, largeobj=30000, minmaxnorm=TRUE) {
  cyto_norm=norm_ch(image=x, index=index, minmaxnorm=minmaxnorm)
  cyto_smooth<-EBImage::filter2(cyto_norm, EBImage::makeBrush(size_smooth, shape="disc"), boundary=c("circular", "replicate"))
  cmask<- mask_ch(cyto_smooth, filter_size, offset, opensize, shape="disc")
  c_sum<-sum_img(cyto_smooth)
  outmask<-EBImage::opening(cyto_smooth>c_sum[2])
  combine <- cmask
  combine[outmask > cmask]=outmask[outmask>cmask]
  combine[y>combine]=y[y>combine]
  cseg<-EBImage::propagate(cyto_smooth, y, lambda=1e-2, mask=combine)
  cseg=EBImage::fillHull(cseg)
  xy<-EBImage::computeFeatures.moment(cseg)[,c('m.cx', 'm.cy')]
  chkpt(xy)
  cf<-EBImage::computeFeatures.shape(cseg)
  cf_Area<-data.frame(cf[,1])
  cf_Area$num=row.names(cf_Area)
  ci=which(cf[,1]>largeobj)
  cseg2=EBImage::rmObjects(cseg, ci, reenumerate=FALSE)

  #! what is this doing?
  xy.nseg<-as.numeric(row.names(EBImage::computeFeatures.moment(y)[,c('m.cx', 'm.cy')]))
  xy.cseg<- as.numeric(row.names(EBImage::computeFeatures.moment(cseg2)[,c('m.cx', 'm.cy')]))
  ind.diff <- setdiff(xy.nseg, xy.cseg)
  nseg<-EBImage::rmObjects(y, ind.diff, reenumerate=F)
  #nseg=reenumerate(nseg) #! I can't find this function do you mean enumerate?
  #cseg=reenumerate(cseg)
  xy.nseg_table<-EBImage::computeFeatures.moment(nseg)[,c('m.cx', 'm.cy')]
  xy.cseg_table<-EBImage::computeFeatures.moment(cseg)[,c('m.cx', 'm.cy')]
  return(list(nuc_seg=nseg, xy_nuc=xy.nseg_table, seg=cseg2, xy_cseg=xy.cseg_table, norm=cyto_norm))
}

###These are Internal Functions invisible to the user and referenced in above (exported) functions###
#' Return single channel from loaded image
#' @inheritParams load_image
#' @return single channel image
norm_ch <- function(image, index, minmaxnorm=TRUE) {
  ch <- image[,,index]
  if (minmaxnorm) {
    minCH = min(as.vector(ch))
    maxCH = max(as.vector(ch))
    ch = EBImage::normalize(ch, ft=c(0,1), c(minCH, maxCH))
  }
  return(ch)
}

#'Make pre-segmented mask for single-channel image
#'@param img processed or normalized single-channel image
#'@param shape shape of brush used to open mask (string)
#'@inheritParams segmentNucleus
mask_ch <- function(img, filter_size, offset, opensize, shape) {
  mask = EBImage::thresh(img, filter_size, filter_size, offset)
  mask = EBImage::opening(mask, EBImage::makeBrush(opensize, shape=shape))
  mask = EBImage::fillHull(mask)
  return(mask)
}

#'Summary of channel intensity, used for image thresholding 
sum_img <- function(img) {
  ar=as.vector(img)
  ar.sum=as.numeric(summary(ar))
  return(ar.sum)
}

#' Checkpoint-warnings
chkpt <- function(maskimg) {
  obj = EBImage::computeFeatures.shape(maskimg)
  if (length(obj)<3) {
    stop("No segmented object was detected, change parameters!")
  }
  if (is.null(nrow(obj))) {
    stop("No segmented object was detected, change parameters!")
  }
  if (nrow(obj)>2000) {
    stop("Poor segmentation, change parameters")
  }
}
chkpt2 <- function(df){
  if (length(xy) < 3 ){
    stop("No segmented objects detected, change setup parameters")
  }
  if (is.null(nrow(xy))){
    stop("No segmented objects detected, change setup parameters")
  }
  if (nrow(xy) > 500 ){
    stop("Poor segmentation, change setup parameters")
  }
}

