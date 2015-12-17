#'A CommonGenesCorPlot function
#'
#' this is a function that creates 4 correlation plots (one per sample)
#' Library "HTqPCR" needed.
#' @param
#' object = qPCR object;
#' @param
#' seqGeneSum = pre processed RNAseq data, for each sample has mean of the log2(CPM).
#' @param
#' extrVal = logical indicator , if TRUE - plots are given for genes with absent signal;
#'                              if FALSE (default value) - plots are given for all the genes in the dataset.
#'
CommonGenesCorPlot<-function(object, seqGeneSum, extrVal=FALSE){
  featureNames(object)<-gsub("_[0-9]*","", featureNames(object))
  rep<-as.numeric(as.factor(gsub(":.+","",colnames(object))))
  geneSum<-t(apply(exprs(object),1,by,rep,mean))
  nSam<-max(rep)
  GenesMatch<-matrix(NA, ncol=nSam*2, nrow=dim(geneSum)[1])
  RNA<-NULL
  for (i in 1:nrow(GenesMatch)) {
    RNA<-seqGeneSum[which(ILM_refseq_gene_COH$Symbol==rownames(geneSum)[i]),]
    GenesMatch[i,]<-c(geneSum[i,], RNA)
  }
  #add row names
  rownames(geneSum)->rownames(GenesMatch)
  #add col names
  name<-NULL
  for (i in 1:nSam) {
    name<-c(name, paste("Taqman",i))
  }
  for (i in 1:nSam) {
    name<-c(name, paste("RNAseq",i))
  }
  colnames(GenesMatch)<-name
  #identify genes that are not in RNAseq data
  notIN<-which(GenesMatch[,1]==GenesMatch[,(nSam+1)])
  GMData<-GenesMatch[-notIN,]
  rownames(GMData)<-rownames(GenesMatch[-notIN,])
  if (extrVal!=TRUE) {
    xlm<-c(round(min(GenesMatch[-notIN,1:nSam]-1),0),35)
    ylm<-c(0,round(max(GenesMatch[-notIN,((1+nSam):(nSam*2))]+0.5),0))
    for (i in 1:nSam){
      plot(GenesMatch[-notIN,i],GenesMatch[-notIN,i+nSam],xlim=xlm, ylim=ylm,
           xlab=paste("TaqMan",i), ylab=paste("RNAseq",i), pch=16, cex=0.6, col=i)
      abline(h=0)
      abline(v=35)
      corp<-cor(GenesMatch[-notIN,i],GenesMatch[-notIN,i+nSam], method="p")
      cors<-cor(GenesMatch[-notIN,i],GenesMatch[-notIN,i+nSam], method="s")
      legend("topright", cex=0.75, legend=c(paste("Pears. cor =",round(corp,3)),
                                            paste("Spearm. cor =",round(cors,3))))
    }
  } else {
    # RNA seq < 0.005
    ind<-which((GMData[,(nSam+1):(nSam*2)]<0.01 | GMData[,1:nSam]>=34.85), arr.ind=TRUE)
    GMData<-GMData[ind,]
    xlm<-c(round(min(GMData[-notIN,1:nSam]-1),0),35)
    ylm<-c(0,round(max(GMData[-notIN,((1+nSam):(nSam*2))]+0.5),0))
    for (i in 1:nSam){
      plot(GMData[-notIN,i],GMData[-notIN,i+nSam],xlim=xlm, ylim=ylm,
           xlab=paste("TaqMan",i), ylab=paste("RNAseq",i), pch=16, cex=0.6, col=i)
      abline(h=0)
      abline(v=35)
      #      corp<-cor(GMData[-notIN,i],GMData[-notIN,i+nSam], method="p")
      #      cors<-cor(GenesMatch[-notIN,i],GenesMatch[-notIN,i+nSam], method="s")
      #      legend("bottomleft", cex=0.75, legend=c(paste("Pears. cor =",round(corp,3)),
      #                                              paste("Spearm. cor =",round(cors,3))))
    }
  }
  return(GMData)
}
