rescale_values_factor_advanced = function(CorrDRSn) {
  
  immse = function(x,y) {
    return(mean((x - y)^2))
  }
  
  fFeat_Corr_v1 = function(x,SMeanFeat1,SMeanFeat2,YcorrDRSn,Crd,i,Feat1idx,Feat2idx) {
    cDRSn=YcorrDRSn
    cDRSn[Crd,i]=YcorrDRSn[Crd,i]*x
    
    CrdFeat1=cDRSn[,3]==Feat1idx
    cMeanFeat1=mean(cDRSn[CrdFeat1,i])
    
    CrdFeat2=cDRSn[,4]==Feat2idx
    cMeanFeat2=mean(cDRSn[CrdFeat2,i])
    
    Res=0.85*immse(cMeanFeat1,SMeanFeat1[i-5])+0.15*immse(cMeanFeat2,SMeanFeat2[i-5]); #0.85 i 0.15 to wspolczynniki dobierane doswiadczalnie
    
    return(Res)
  }
  
  
  YcorrDRSn=CorrDRSn
  MaxC=0
  
  for i in 1:size(CorrDRSn,2)-5 { #bo dane od 6 kolumny
    for j in 1:max(YcorrDRSn[,3]) { #Czyli po cesze 1
      CrdFeat1=YcorrDRSn[,3]==j
      if (length(CrdFeat1)==0)
        YcorrMeanFeat1[i,j]=0
      else
        YcorrMeanFeat1[i,j]=mean(YcorrDRSn[CrdFeat1,i+5])
    }

    for j in 0:max(YcorrDRSn[,4]) { #czyli po cesze 2
      CrdFeat2=YcorrDRSn[:,4]==j
      if (length(CrdFeat2)==0) 
        YcorrMeanFeat2[i,j+1]=0
      else
        YcorrMeanFeat2[i,j+1]=mean(YcorrDRSn[CrdFeat2,i+5])
    }
  }
  SMeanFeat1=mean(YcorrMeanFeat1[:,3:22])
  SMeanFeat2=mean(YcorrMeanFeat2)
  
  for i in 6:ncol(CorrDRSn) {         #Dane od 6 kolumny!
    for Feat1idx=1:max(YcorrDRSn[,3]) {
      for Feat2idx=0:max(YcorrDRSn[,4]) {
        Crd=YcorrDRSn[:,3]==Feat1idx & YcorrDRSn[:,4]==Feat2idx
        if (length(Crd)>0) 
          CorrCoeff = fminbnd(@(x) fFeat_Corr_v1(x,SMeanFeat1,SMeanFeat2,YcorrDRSn,Crd,i,Feat1idx,Feat2idx),0,1000000)
        if CorrCoeff>MaxC
          MaxC=CorrCoeff
        YcorrDRSn[Crd,i]=YcorrDRSn[Crd,i]*CorrCoeff
      }
    }
  }

  return(YcorrDRSn)
}
