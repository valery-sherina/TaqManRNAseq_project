#'A qPCRdata function
#'
#' this is a function that makes a new qPCR object. Library "HTqPCR" needed.
#' @param
#' DataName = a string with the name of the data;
#' @param
#' ControlGene = the name of the control gene.
#'
qPCRdata<-function(DataName, ControlGene) {
  object<-get(DataName)
  #  ind<-as.numeric(gsub("\\w+_","",DataName))

  ft <- rep("Target",nrow(object))
  ind<-which(rownames(object)==ControlGene)
  ft[ind] <- "Endogenous Control"

  fc <- matrix("OK",nrow=nrow(object),ncol=ncol(object))
  fc[which(object>34.9999,arr.ind=TRUE)] <- "Undetermined"
  colnames(fc) <- colnames(object)
  rownames(fc) <- rownames(object)

  fl <- matrix("Passed",nrow=nrow(object),ncol=ncol(object))
  fl[which(object>34.9999,arr.ind=TRUE)] <- "Flagged"
  colnames(fl) <- colnames(object)
  rownames(fl) <- rownames(object)

  TaqMan <- new("qPCRset", exprs=object, flag=fl)
  featureNames(TaqMan) <- rownames(object)
  featureType(TaqMan) <- ft
  featureCategory(TaqMan) <- as.data.frame(fc)

  sType <- c(rep("A",4),rep("B",4),rep("C",4),rep("D",4))
  #  PartID <- c(rep(ind, dim(object)[2]))
  #  tab <- data.frame(sampleName=colnames(object), sampleType=sType, dataPartID=PartID)
  tab <- data.frame(sampleName=colnames(object), sampleType=sType)
  phenoData(TaqMan) <- AnnotatedDataFrame(data=tab)
  return(TaqMan)
}
