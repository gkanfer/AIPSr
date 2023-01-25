#' Loads Image
#' This function loads an image from file-name, then normalizes each channel
#' @param file_name path to image file (string)
#' @return Image object (EBImage)
#' @export

LoadImage <- function(file_name) {
  x<-EBImage::readImage(file_name)
  return(x)
}
