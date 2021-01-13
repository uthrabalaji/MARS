#rm(list=ls())
#suppressPackageStartupMessages(library(DESeq2))
if("DESeq2" %in% rownames(installed.packages()) == FALSE) {install.packages("DESeq2")} 
suppressPackageStartupMessages(library(DESeq2))

MARS <- function(Rminus_BSJ,Rminus_GC,Rplus_BSJ,design_file)
{
  #Ribominus/RNaseR-
  sf.gc.rminus = estimateSizeFactorsForMatrix(counts=Rminus_GC)
  keep = rowMeans(Rminus_BSJ)>1
  Rminus_BSJ.filt = Rminus_BSJ[keep,]
  dds.rminus <- DESeqDataSetFromMatrix(Rminus_BSJ.filt,colData = design_file,design = ~0+Group)
  sizeFactors(dds.rminus) <- sf.gc.rminus
  dds.rminus <- DESeq(dds.rminus,quiet = T)
  res.rminus <- as.data.frame(results(dds.rminus,alpha = 0.05))
  
  #Ribominus/RNaseR+
  sf.bsj.rplus = estimateSizeFactorsForMatrix(counts = Rplus_BSJ)
  sf.bsj.rminus = estimateSizeFactorsForMatrix(counts = Rminus_BSJ)
  #compound factor
  sf.compound = (sf.bsj.rplus * sf.gc.rminus)/sf.bsj.rminus
  keep = rowMeans(Rplus_BSJ)>1
  Rplus_BSJ.filt = Rplus_BSJ[keep,]
  dds.rplus <- DESeqDataSetFromMatrix(Rplus_BSJ.filt,colData = design_file,design = ~0+Group)
  sizeFactors(dds.rplus) <- sf.compound
  dds.rplus <- DESeq(dds.rplus,quiet = T)
  res.rplus <- as.data.frame(results(dds.rplus,alpha = 0.05))
  
  #meta-analysis 
  common.circ = intersect(rownames(res.rminus),rownames(res.rplus))
  unique.rminus = rownames(res.rminus)[-match(common.circ,rownames(res.rminus))]
  unique.rplus = rownames(res.rplus)[-match(common.circ,rownames(res.rplus))]
  
  res.rminus.common = res.rminus[common.circ,]
  res.rplus.common = res.rplus[common.circ,]
  
  p.rminus = res.rminus.common$pvalue; names(p.rminus) = rownames(res.rminus.common)
  p.rplus = res.rplus.common$pvalue; names(p.rplus) = rownames(res.rplus.common)
  
  # check for NA in both rminus and rplus
  na.rminus = names(p.rminus)[which(is.na(p.rminus))]
  na.rplus = names(p.rplus)[which(is.na(p.rplus))]
  na.common = intersect(na.rminus,na.rplus)
  unique.rminus = c(unique.rminus,na.rplus[-match(na.common,na.rplus)])
  unique.rplus = c(unique.rplus,na.rminus[-match(na.common,na.rminus)])
  
  common.circ = common.circ[-match(c(na.common,na.rminus,na.rplus),common.circ)]
  res.rminus.common = res.rminus[common.circ,]
  res.rplus.common = res.rplus[common.circ,]
  
  #calculate z-score
  p.rminus = res.rminus.common$pvalue; names(p.rminus) = rownames(res.rminus.common)
  p.rplus = res.rplus.common$pvalue; names(p.rplus) = rownames(res.rplus.common)
  
  lfc.rminus = res.rminus.common$log2FoldChange; names(lfc.rminus) = rownames(res.rminus.common)
  lfc.rplus = res.rplus.common$log2FoldChange; names(lfc.rplus) = rownames(res.rplus.common)
  lfc.avg = apply(data.frame(lfc.rminus,lfc.rplus),1,mean)
  
  z.rminus = qnorm(p.rminus/2,lower.tail = F)*sign(lfc.rminus)
  z.rplus = qnorm(p.rplus/2,lower.tail = F)*sign(lfc.rplus)
  
  #weights
  lfcSE.rminus = res.rminus.common$lfcSE; names(lfcSE.rminus) = rownames(res.rminus.common)
  lfcSE.rplus = res.rplus.common$lfcSE; names(lfcSE.rplus) = rownames(res.rplus.common)
  w.rminus = 1/lfcSE.rminus
  w.rplus = 1/lfcSE.rplus
 
  #covariance
  idx1 = which(abs(z.rminus)<1.96)
  idx2 = which(abs(z.rplus)<1.96)
  idx = intersect(idx1,idx2)
  cov.out = cov(z.rminus[idx],z.rplus[idx])

  df = data.frame(z.rminus=z.rminus,z.rplus=z.rplus,w.rminus=w.rminus,w.rplus=w.rplus,w.rplus,cov.out=cov.out)
  z.meta =  apply(df,1,function(x){((x["w.rminus"]*x["z.rminus"])+(x["w.rplus"]*x["z.rplus"]))/sqrt((x["w.rminus"]^2)+(x["w.rplus"]^2)+(2*x["w.rminus"]*x["w.rplus"]*x["cov.out"]))})
  p.meta <- 2 * (1 - pnorm(abs(z.meta)))
  
  p.unique.rminus <- res.rminus$pvalue[match(unique.rminus,rownames(res.rminus))]
  names(p.unique.rminus) <- unique.rminus
  p.unique.rminus <- p.unique.rminus[!is.na(p.unique.rminus)]
  lfc.unique.rminus <- res.rminus$log2FoldChange[match(names(p.unique.rminus),rownames(res.rminus))]
  names(lfc.unique.rminus) <- names(p.unique.rminus)
  
  p.unique.rplus <- res.rplus$pvalue[match(unique.rplus,rownames(res.rplus))]
  names(p.unique.rplus) <- unique.rplus
  p.unique.rplus <- p.unique.rplus[!is.na(p.unique.rplus)]
  lfc.unique.rplus <- res.rplus$log2FoldChange[match(names(p.unique.rplus),rownames(res.rplus))]
  names(lfc.unique.rplus) <- names(p.unique.rplus)
  
  all.circ = c(names(p.meta),names(p.unique.rminus),names(p.unique.rplus))
  res.all.circ = cbind(res.rminus[all.circ,c("baseMean","log2FoldChange","pvalue")],res.rplus[all.circ,c("baseMean","log2FoldChange","pvalue")])
  rownames(res.all.circ) = all.circ
  colnames(res.all.circ) = paste(rep(c("Ribominus/RNaseR-","Ribominus/RNaseR+"),each=3),rep(c("baseMean","log2FoldChange","pvalue"),2))
  
  meta.out <- data.frame(circRNA = all.circ , res.all.circ, avg.log2FC = c(lfc.avg,lfc.unique.rminus,lfc.unique.rplus),meta.pvalue = c(p.meta,p.unique.rminus,p.unique.rplus),check.names = F)
  
  return(meta.out)
}


