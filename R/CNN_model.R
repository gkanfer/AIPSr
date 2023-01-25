## 3 CNN Training

#! Should CNN training have its own pickCells (classification-based selection) function

#' 3.1 Generation and upload of single cell images for CNN training
#' @param 
#' @param pheno0 String. Name of control (normal) phenotype
#' @param pheno1 String. Name of positive phenotype
#' @param label_class String. Classification of negative or positive phenotype?
#' @return 
#' @export
genImageSet <- function(cseg, base_dir, pheno0, pheno1, int){
  dir.create(file.path(base_dir, pheno0))
  dir.create(file.path(base_dir, pheno1))
  nseg = cseg$nuc_seg
  x=cseg$norm
  #classify cells 
  seg_disp = EBImage::paintObjects(nseg, EBImage::toRGB(x*int),opac=c(1,0.8),col=c("Green",NA), thick=TRUE,closed=FALSE)
  EBImage::display(seg_disp,"raster")
  celltext = text(x=xy_nuc[,1], y=xy_nuc[,2], labels="", col="yellow", cex=font_size)
  c<-0
  readline(paste0("Select Positive cells"))
  temp<-locator()
  c<-c(c,temp)
  xy<-EBImage::computeFeatures.moment(nseg)[,c('m.cx','m.cy')]
  df.c <- cbind(c$x,c$y)
  knn.out <- yaImpute::ann(as.matrix(xy), as.matrix(df.c), k=2)
  pos_row<-knn.out$knnIndexDist
  
  readline(paste0("Select Negative cells"))
  temp<-locator()
  c<-c(c,temp)
  xy<-EBImage::computeFeatures.moment(nseg)[,c('m.cx','m.cy')]
  df.c <- cbind(c$x,c$y)
  knn.out <- yaImpute::ann(as.matrix(xy), as.matrix(df.c), k=2)
  neg_row<-knn.out$knnIndexDist
  
  #Save individual cells/files in phenotype dir
  stack<-EBImage::stackObjects(cseg$seg, x)
  
  if (pos_row != 0) {
    p_frames=getFrames(stack, i=pos_row)
    for (i in 1:length(p_frames)) {
      img = p_frames[[i]]
    }
  }
  if (neg_row !=0) {
    n_frames=getFrames(stack, i=neg_row)
    for (i in 1:length(n_frames)) {
      
    }
  }
  
  
  #Break the 20x magnification image into single cell images
  #stack<-EBImage::stackObjects(cseg$seg, x,ext = c(satck_size, satck_size)) 
  #! do we need the ext option? Where does satck_size come from?
  
  #! I don't think we should include this preprocessing function. Move to generator creation if we still want that option
  # Fimg <- function(x) {
  #   x1<-readImage(x, type="tiff")
  #   dx <- (dim(x1)[1])
  #   dx1 <- dx%/%8
  #   ft <- matrix(c(-0.05,-0.02,0, -0.02, -0.05), nrow=5, ncol=5)
  #   ft[5,5] <- 2
  #   y1 <- filter2(x1, ft)
  #   y2 <- y1[(dx1-10):(dx-dx1),(dx1-10):(dx-dx1)]
  #   writeImage(y2, paste(x),type = "tiff", bits.per.sample = 8L, compression = "LZW")
  # }
  
}

#' 3.2 Create testing, training directories (split based off of user-inputed ratio) and move files to appropriate destinations
#' Create files into testing and training sets, create test directory. Remaining files in original base_dir=train_dir 
#' @param base_dir 
#' @param pheno0 String. Name of control (normal) phenotype
#' @param pheno1 String. Name of positive phenotype
#' @param test_ratio Number between 0-1. Proportion of images to relocate to testing directory
#' @export
splitDirs <- function(base_dir, pheno0, pheno1, test_ratio) {
  cls_list = c(pheno0, pheno1)
  for (cls in cls_list) {
    createdir(file.path(base_dir,'test',cls))
    createdir(file.path(base_dir,'train',cls))
    src = file.path(base_dir, cls) 
    allfiles = list.files(path=src, full.names=TRUE)
    allfiles = sample(allfiles)  #random shuffling 
    test_n= as.integer(length(allfiles)*test_ratio)
    test_files=allfiles[0:test_n]
    train_files=allfiles[test_n:as.integer(length(allfiles))]
    print(paste0("Testing-",cls,": ",length(test_files)))
    print(paste0("Training-",cls,": ",length(train_files)))
  }
  for (i in test_files) {
    path = file.path(base_dir,"test",cls) 
    file.copy(i, path)
  }
}

##Internal Functions##
#'@param  cseg Output (list) from segmentCyto
#'


#' #' Directory create
#' createdir <- function(path, subpath) {
#'   dir <- file.path(path, subpath)
#'   dir.create(dir)
#'   return(dir)
#' }
