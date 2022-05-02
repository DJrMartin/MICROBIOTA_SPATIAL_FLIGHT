#####FONCTION POUR LA NORMALISATION
library(metagenomeSeq)
NORMALISATION_MICROBIOTA <- function(data, approches="BRUTES"){
  OTU_normalised=NULL
  if (approches=="BRUTES"){
    OTU_normalised=as.data.frame(data)
  }
  if (approches=="TSS"){#relative abundances
    for (ind in 1:dim(data)[1]){
      ID=data[ind,]/sum(data[ind,])
      OTU_normalised=rbind(OTU_normalised,ID)
    }
  }
  if (approches=="CLR"){#Center log ratio
    for (ind in 1:dim(data)[1]){
      x=data[ind,which(matrix_otu[ind,]>0)]
      ID=log(x/(prod(x)^(1 / length(x))))
      OTU_normalised=rbind(OTU_normalised,ID)
    }
  }
  if (approches=="CSS"){#median-like quantile normalisation
    ##Paulson, J.N., Stine, O.C., Bravo, H.C. & Pop, M. Nat. Methods 10, 1200–1202 (2013)
    metaSeqObject=newMRexperiment(as.data.frame(t(data)))
    metaSeqObject_CSS = cumNorm(metaSeqObject, 
                                p=cumNormStatFast(metaSeqObject))
    OTU_normalised=data.frame(t(MRcounts(metaSeqObject_CSS, norm=TRUE, log=TRUE)))
  }
  OTU_normalised=data.frame(OTU_normalised)
  return(OTU_normalised)
}

  