#' 2.1 Measures cell features
#' @param x (Normalized) image of cytosolic/phenotype channel. Value returned from segmentCyto (out$norm) or user generated via EBImage. 
#' @param mask_cyto Cell segmentation mask output from segmentCyto
#' @param label_class String. Positive or Negative. 
#' @return List of table, and row-names table (?) with features (cytosolic img) for each segmented cell 
#' @export

extractFeatures <- function(x, mask_cyto, label_class="Positive") {
  table_shape = EBImage::computeFeatures.shape(mask_cyto,x)
  table_moment = EBImage::computeFeatures.moment(mask_cyto,x)
  table_basic = EBImage::computeFeatures.basic(mask_cyto,x)
  table_test <- as.data.frame(cbind(table_basic,table_moment,table_shape))
  #table_test$predict <- label_class 
  ##! what is purpose of rownames, table restructuring ?
  rownameTable<-row.names(table_test)
  #table_test$predict<-label_class
  feature_table<-data.frame(cbind(rownameTable,table_test))
  Ts.mix<-feature_table[,2:20]
  rowNameTable<-feature_table[,1]
  Ts.mix$predict<-label_class #! what is the point of this? Won't the predict column be filled in during pickCells function?
  Features_Table<-list()
  Features_Table[["Ts.mix"]]<-Ts.mix
  Features_Table[["rowNameTable"]]<-rowNameTable
  return(Features_Table)
}

#' 2.2 Choose cells for classification
#' @param
#' @return
#' @export
pickCells<-function(mask_nuc, x, xy_nuc, Ts.mix, int, font_size=0.7, label_class="Postive", display_select=TRUE) {
  seg_disp = EBImage::paintObjects(mask_nuc, EBImage::toRGB(x*int),opac=c(1,0.8),col=c("Green",NA),thick=TRUE,closed=FALSE)
  EBImage::display(seg_disp,"raster")
  celltext = text(x=xy_nuc[,1], y=xy_nuc[,2], labels="", col="yellow", cex=font_size)
  c<-0
  readline(paste0("Select", label_class, "cells"))
  temp<-locator()
  c<-c(c,temp)
  xy<-EBImage::computeFeatures.moment(mask_nuc)[,c('m.cx','m.cy')]
  df.c <- cbind(c$x,c$y)
  knn.out <- yaImpute::ann(as.matrix(xy), as.matrix(df.c), k=2)
  row_numb<-knn.out$knnIndexDist
  class(row_numb)
  row_numb<-as.data.frame(row_numb)
  Ts.mix$predict<-0
  if (label_class=="Postive"){
    Ts.mix[row_numb$V1,20]<-"P"
    Ts.mix[-row_numb$V1,20]<-"N"
  }else{
    Ts.mix[row_numb$V1,20]<-"N"  
    Ts.mix[-row_numb$V1,20]<-"P"
  }
  nr = which(Ts.mix$predict %in% "N")
  
  Table_class<-list()
  if (display_select) {
    seg.temp = EBImage::rmObjects(mask_nuc, nr)
    seg_class_sel = EBImage::paintObjects(seg.temp,EBImage::toRGB(x*int),opac=c(1,0.8),col=c("Red",NA),thick=TRUE,closed=FALSE)
    EBImage::display(out$segsel)
  }
  return(list(table_train=Ts.mix, segsel=seg_class_sel))
} 

##! How to train with multiple images? Let user make their own loop using segmentNuc -> segmentCyto -> pickCells ? Loop would need to return a single table for all images in order to work with current modelSVM function

##! Does SVM model function incoroporate both positive and negative classification? Should there be two table arguments (pos, neg), or should their pickCells output actually be a single table for both negative and positive.. but how does that work? 

#'2.3 SVM model 
#' @param total_train Combined table of both positive and negative classified cells
#' @param kernel_linear ??
#' @param cost ?
#' @param degree ?
#' @return model accuracy **and model??
#' @export
modelSVM<-function(total_train, kernel_linear=TRUE,cost= 10, degree = 45){
  ind0<-which(is.na(total_train$predict))
  if (length(ind0)>1) total_train<-total_train[-ind0,]
  ind<-which(is.na(total_train$predict))
  if (length(ind)>1) total_train<-total_train[-ind,]
  ind1<-grep("\\d",total_train$predict)
  if (length(ind)>1) total_train<-total_train[-ind1,]  
  x<-total_train[,2:19]
  y<-total_train[,20]
  if (kernel_linear){
    acc<-rep(0,100)
    pb <- progress::progress_bar$new(total = 100)
    for (i in 1:length(acc)){
      pb$tick()
      TestIndex<-sample(nrow(x),round(nrow(x)/2))
      model<-e1071::svm(x=x[TestIndex,], y=as.factor(y[TestIndex]),kernel="linear",type = "C",cost= 10, degree = 45,probability = TRUE)
      y.pred<-e1071::predict(model,x[-TestIndex,])
      acc[i]=length(which(y.pred==y[-TestIndex]))/length(y.pred)
    }
  }else{
    acc<-rep(0,100)
    pb <- progress::progress_bar$new(total = 100)
    for (i in 1:length(acc)){
      pb$tick()
      TestIndex<-sample(nrow(x),round(nrow(x)/2))
      model<-e1071::svm(x=x[TestIndex,], y=as.factor(y[TestIndex]),kernel="radial",type = "C",cost= 10, degree = 45,probability = TRUE)
      y.pred<-e1071::predict(model,x[-TestIndex,])
      acc[i]=length(which(y.pred==y[-TestIndex]))/length(y.pred)}
  }
  mean(acc)
  #return(acc) #! should this function also return a model?
  return(list(model=model, acc=acc))
}

#' 2.4 PCA Analysis: features reduction 
#' @param 
#' @return 
#' @export

featuresPCA<-function(Table_class_train, cartesian_lim_x=c(-12,10), cartesian_lim_y=c(-10,5), font_size=14){
  ind0<-which(is.na(Table_class_train$predict))
  if (length(ind0)>1) Table_class_train<-Table_class_train[-ind0,]
  ind<-which(is.na(Table_class_train$predict))
  if (length(ind)>1) Table_class_train<-Table_class_train[-ind,]
  ind1<-grep("\\d",Table_class_train$predict)
  if (length(ind)>1) Table_class_train<-Table_class_train[-ind1,]   
  Table_class_train<-Table_class_train %>% filter_at(vars(1:20), all_vars(!is.infinite(.)))
  x<-1
  z<-sample(1:nrow(Table_class_train),nrow(Table_class_train))
  while (x<10){
    myPr <- try(prcomp(Table_class_train[z, 2:19], scale = TRUE))
    if (inherits(myPr, "try-error")) {
      x<-x+1
      z<-sample(1:nrow(Table_class_train),round(nrow(Table_class_train)/x,2))
    }else{
      break
    }
  }
  PCA_PLOT<-ggbiplot(myPr, pc.biplot = T,obs.scale = 2, var.scale = 1, groups = F, ellipse = F, circle = F,alpha = 0.02,varname.adjust = 5)+
    coord_cartesian(xlim = cartesian_lim_x,ylim = cartesian_lim_y)+
    theme(axis.title.x =element_text(size = font_size))+
    theme(axis.title.y =element_text(size = font_size))+
    theme(axis.text = element_text(size=font_size))+
    theme(axis.line = element_line(colour = "black",size=1), panel.border = element_blank(),panel.background = element_blank())+
    labs(title="PCA - feature selection") + 
    theme(plot.title = element_text(hjust = 0.5))+
    theme(panel.border = element_blank(),panel.background = element_blank())+
    theme(panel.grid.minor = element_line(colour = "white"))+
    theme(panel.border = element_blank(),panel.background = element_blank())+
    theme(panel.grid.minor = element_line(colour = "white"))+
    theme(legend.position="none")
  return(PCA_PLOT)
}

#2.4 Features selection and SVM model tuning
#' @param 
#' @return 
#' @export
SVM_FeatureSel<-function(Table_class_train,kernel_linear=TRUE,cost= 10, degree = 45){
  ind0<-which(is.na(Table_class_train$predict))
  if (length(ind0)>1) Table_class_train<-Table_class_train[-ind0,]
  ind<-which(is.na(Table_class_train$predict))
  if (length(ind)>1) Table_class_train<-Table_class_train[-ind,]
  ind1<-grep("\\d",Table_class_train$predict)
  if (length(ind)>1) Table_class_train<-Table_class_train[-ind1,]   
  Table_class_train<-Table_class_train %>% filter_at(vars(1:20), all_vars(!is.infinite(.)))
  vec<-NULL  
  for (i in 2:19){
    print(colnames(Table_class_train[i]))
    vec.temp<-readline("Include Feture ? Y/N")
    vec.temp<-as.character(vec.temp)
    if (vec.temp=="Y"){
      vec<-c(vec,i)
    }
    rm(vec.temp)  
  }
  Table_class_sel<-Table_class_train[,c(vec,which(colnames(Table_class_train)=="predict"))]
  x<-Table_class_sel[,-(ncol(Table_class_sel))]
  y<-Table_class_train$predict
  pb <- progress_bar$new(total = 100)
  if (kernel_linear){
    acc<-rep(0,100)
    for (i in 1:length(acc)){
      pb$tick()
      TestIndex<-sample(nrow(x),round(nrow(x)/2))
      model<-svm(x=x[TestIndex,], y=as.factor(y[TestIndex]),kernel="linear",type = "C",cost= 10, degree = 45,probability = TRUE)
      y.pred<-predict(model,x[-TestIndex,])
      acc[i]=length(which(y.pred==y[-TestIndex]))/length(y.pred)
    }
  }else{
    acc<-rep(0,100)
    for (i in 1:length(acc)){
      pb$tick()
      TestIndex<-sample(nrow(x),round(nrow(x)/2))
      model<-svm(x=x[TestIndex,], y=as.factor(y[TestIndex]),kernel="radial",type = "C",cost= 10, degree = 45,probability = TRUE)
      y.pred<-predict(model,x[-TestIndex,])
      acc[i]=length(which(y.pred==y[-TestIndex]))/length(y.pred)}
  }
  model_svm<-list()
  model_svm[["ACC"]]<-acc
  model_svm[["Selected_Features"]]<-colnames(x)
  model_svm[["model"]]<-model
  return(model_svm)
}
