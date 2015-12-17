#'A toGetData function
#'
#' this is a function that subsets the data and back-calculates
#' Ct values, it returns vector of Ct values
#' @param
#' data = vector of qPCR signal;
#' @param
#' NDstatus = vector of nondetect status;
#' @param
#' trtSubs = the first position of different treatments,
#' here it is 1,5,9 or 13;
#' @param
#' from = the first record og the plate,
#' here it is 1, 97, 193, 289, 385, 481, 577, 673, 769, 865, 961
toGetData<-function(data, NDstatus, from, trtSubs){
  to<-from+95
  if (to > dim(data)[1]) { to<-dim(data)[1] }
  #sunset the data
  temp <-data[from:to,  trtSubs:(trtSubs+3)]
  #get Ct for control
  Ct_control<-log2(unique(data[from:to,  trtSubs:(trtSubs+3)][NDstatus[from:to,  trtSubs:(trtSubs+3)]==1]))+35
  #calculate Ct values
  Ct_values<- Ct_control - log2(temp)
  return(Ct_values)
}
