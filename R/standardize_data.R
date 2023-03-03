#' Read and standardize microarray data
#'
#' @param targets data frame with experiment information, obtained using limma::readTargets function
#' @param correct conduct additional factor based data correction (TRUE/FALSE)
#'
#' @return
#' Normalized expression table  
#'
#' @export
#' @examples
#' targets <- readTargets("sample_desc.txt",path=path)
#' exprs_correctted = standardize_data(targets,TRUE)
standardize_data = function(targets,correct=FALSE) {
  
  RG <- read.maimages(targets, columns = list(G = "gMedianSignal", Gb = "gBGMedianSignal", R = "rMedianSignal",
                                              Rb = "rBGMedianSignal"), annotation = c("Row", "Col","FeatureNum", "ControlType","ProbeName","SystematicName"),path=path)
  
  RG.bg <- backgroundCorrect(RG, method="minimum", offset=1)
  
  ### correction ###
  if (correct) {
    tFval = RG.bg$genes$ProbeName
    tFval_mat = strsplit2(tFval,'_')
    factor_values = tFval_mat[,4]
    factor_values[tFval_mat[,2]!="GC"] = NA
    factor_values = as.numeric(factor_values)
    RG.bg$G = rescale_values_factor_simple(RG.bg$G,factor_values)
    RG.bg$R = rescale_values_factor_simple(RG.bg$R,factor_values)
  }
  
  MA <- normalizeWithinArrays(RG.bg, method="loess")
  MA <- normalizeBetweenArrays(MA, method="Aquantile")
  MA.avg <- avereps(MA, ID=MA$genes$SystematicName)
  
  RG.pq <- RG.MA(MA)
  exptab <- cbind(RG.pq$R, RG.pq$G)
  colnames(exptab)=c(paste(MA$targets$Cy3,'R',sep='_'), paste(MA$targets$Cy5,'G',sep='_'))
  exptabAnnot = cbind(MA$genes[,c('ProbeName','SystematicName')],exptab)
  
  ids = strsplit2(exptabAnnot$ProbeName,'_')
  exptabAnnot$Feature = ids[,4]
  exptabAnnot$ProbeID = paste(ids[,1],ids[,2],ids[,3],sep="_")
  exptabAnnot$ProbeType = ids[,2]
  exptabAnnotSub = exptabAnnot[ids[,2] %in% c('GC','LN','An','G4','T7'),]
  
  exptabAnnotSub = exptabAnnotSub[order(exptabAnnotSub$ProbeID),]
  SampleIDs = colnames(exptab)
  
  annotCols = colnames(exptabAnnotSub)[!colnames(exptabAnnotSub) %in% SampleIDs]
  annotCols = c('ProbeName','SystematicName','ProbeID','ProbeType','Feature')
  exptabAnnotSub_export = exptabAnnotSub[,c(annotCols,SampleIDs)]
  
  
  CombinedProbeData = readRDS('../pProbes_CombinedData_v1.RDS')
  AnnotColumns = mcols(CombinedProbeData[,c('Probe.nID','Exon.exon_id')])
  
  exptabAnnotSub_export2 = merge(AnnotColumns,exptabAnnotSub_export,by.x='Probe.nID',by.y='ProbeName')
  
  return(exptabAnnotSub_export2)
}