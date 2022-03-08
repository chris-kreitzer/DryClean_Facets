## Traditional Facets algorithm;
## minor modifications to run the quality assessment
## 
## 03/07/2022
## chris-kreitzer
##

facets_original = function(count_matrix, cval){
  require(facets)
  require(pctGCdata)
  source('~/Documents/GitHub/DryClean_Facets/Scripts/FacetsQC.R')
  count_matrix = facets::readSnpMatrix(count_matrix)
  data_proc = facets::preProcSample(rcmat = count_matrix, ndepth = 35, ndepthmax = 1000, gbuild = 'hg19')
  out = facets::procSample(x = data_proc, cval = cval)
  fit = facets::emcncf(out)
  
  fit$cncf = cbind(fit$cncf, cf = out$out$cf, tcn = out$out$tcn, lcn = out$out$lcn)
  fit$cncf$lcn[fit$cncf$tcn == 1] = 0
  fit$cncf$lcn.em[fit$cncf$tcn.em == 1] = 0
  
  # Generate output
  full_out = list(snps = out$jointseg,
                  segs = fit$cncf,
                  purity = as.numeric(fit$purity),
                  ploidy = as.numeric(fit$ploidy),
                  dipLogR = out$dipLogR,
                  alBalLogR = out$alBalLogR,
                  flags = out$flags,
                  em_flags = fit$emflags,
                  loglik = fit$loglik)
  
  qc = facets_fit_qc(facets_output = full_out)
  
  return(list(facets_output = full_out,
              qc = qc))
}
  
# y = facets_original(count_matrix = '~/Desktop/countsMerged____P-0001396-T05-IM6_P-0001396-N01-IM6.dat.gz', cval = 150)
# y$qc$facets_qc
