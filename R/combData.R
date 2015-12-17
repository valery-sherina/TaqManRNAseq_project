#' A combData function
#'
#' this is a function that makes a dataset for one set
#' of genes (one plate)
#' @param
#' deltaCt = delta Ct values;
#' @param
#' ND = vector of nondetect status;
#' @param
#' sbsBgn = the first position of different treatments;
#' @param
#' gene = gene names;
#' @param
#' delNondet = an indicator whether or not we delete genes
#' that are unexpressed in any of the samples, default = TRUE

combData<-function(deltaCt,ND,sbsBgn,gene, delNondet=TRUE){
  # subset the data and back-calculate Ct values
  tst1<-toGetData(data=deltaCt,NDstatus=ND,from=sbsBgn,trtSubs=1)
  tst5<-toGetData(deltaCt,ND,sbsBgn,5)
  tst9<-toGetData(deltaCt,ND,sbsBgn,9)
  tst13<-toGetData(deltaCt,ND,sbsBgn,13)

  # combine in one data set with 4 samples
  tst<-cbind(tst1,tst5,tst9,tst13)

  # add names for rows and columns
  sbsEnd<-sbsBgn+95
  if ( sbsEnd > dim(deltaCt)[1] ) { sbsEnd<-dim(deltaCt)[1] }
  rownames(tst)<-gene[sbsBgn:sbsEnd]
  colnames(tst)<-c("samp1:1","samp1:2","samp1:3","samp1:4",
                   "samp2:1","samp2:2","samp2:3","samp2:4",
                   "samp3:1","samp3:2","samp3:3","samp3:4",
                   "samp4:1","samp4:2","samp4:3","samp4:4")
  if (delNondet == TRUE) {
    # delete nondetected genes
    allND<-c(which(apply(tst1,1,sum)==140), which(apply(tst5,1,sum)==140),
             which(apply(tst9,1,sum)==140), which(apply(tst13,1,sum)==140))
    allND<-unique(sort(allND))
    tst<-tst[-allND,]
  }
  return(tst)
}
