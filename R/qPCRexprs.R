#'A qPCRexprs function
#'
#' this is a function that makes one qPCR class object from two qPCR objects of one kind. Library "HTqPCR" needed.
#' @param
#' object1 = first object;
#' @param
#' object2 = second object.
#'
qPCRexprs <- function(object1, object2) {
  ft <- c(featureType(object1),featureType(object2))

  fc <- rbind(featureCategory(object1),featureCategory(object2))

  fl <- matrix("Passed",nrow=nrow(fc),ncol=ncol(fc))
  fl[which(fc!="OK",arr.ind=TRUE)] <- "Flagged"
  colnames(fl) <- colnames(fc)
  rownames(fl) <- rownames(fc)

  gexprs <- rbind(exprs(object1),exprs(object2))
  TaqMan <- new("qPCRset", exprs=gexprs, flag=fl)
  featureNames(TaqMan) <- c(featureNames(object1),featureNames(object2))
  featureType(TaqMan) <- ft
  featureCategory(TaqMan) <- as.data.frame(fc)

  sType <- c(rep("A",4),rep("B",4),rep("C",4),rep("D",4))
  #  PartID <- c(rep(ind1, dim(object1)[2]),rep(ind2, dim(object2)[2]))
  tab <- data.frame(sampleName=pData(object2)$sampleName,
                    sampleType=pData(object2)$sampleType)
  phenoData(TaqMan) <- AnnotatedDataFrame(data=tab)
  return(TaqMan)
}
