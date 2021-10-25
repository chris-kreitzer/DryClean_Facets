## Facets compilation file:
## This file conditionally replaces CnLR estimates from Facets with those provided by DryClean
## Conditionally, because not all positions are available from DryClean

run_facets_cleaned = function(read_counts,
                              read_cleaned,
                              MODE = c('full', 'partial'),
                              cval = 100,
                              dipLogR = NULL,
                              ndepth = 35,
                              snp_nbhd = 250,
                              min_nhet = 15,
                              seed = 100) {
  message('MODE must either be partial, or full')
  
  if(MODE == 'partial'){
    print(paste0('Mode: ', MODE, ' is selected'))
    
    # Check input 
    missing_cols = setdiff(c('Chromosome', 'Position', 'NOR.DP', 'TUM.DP', 'NOR.RD', 'TUM.RD'), names(read_counts)) 
    if (length(missing_cols) > 0) {
      stop(paste0('Input missing column(s)', paste(missing_cols, collapse = ', '), '.'), call. = FALSE)
    }
    
    set.seed(seed = seed)
    genome = 'hg19'
    
    # Run FACETS algorithm
    dat = facets::preProcSample(rcmat = read_counts, 
                                ndepth = ndepth, 
                                ndepthmax = 1000,
                                het.thresh = 0.25, 
                                snp.nbhd = snp_nbhd, 
                                cval = 25,
                                gbuild = genome, 
                                hetscale = TRUE, 
                                unmatched = FALSE)
    
    
    #' substitute cnLR from original Facets run with DryClean input
    #' DryClean
    data_cleaned = read_cleaned
    data_cleaned$bin = paste(data_cleaned$seqnames, data_cleaned$start, sep = ';')
    #' Facets
    preProc_jointseg = dat$jointseg
    preProc_jointseg$bin = paste(preProc_jointseg$chrom, preProc_jointseg$maploc, sep = ';')
    missing_cnlr = which(is.na(preProc_jointseg$cnlr))
    
    #' replace original CnLR from Facets with DryClean's foreground.log (where applicable)
    ii = which(preProc_jointseg$bin %in% data_cleaned$bin)
    jj = which(data_cleaned$bin %in% preProc_jointseg$bin[ii])
    
    preProc_jointseg$cnlr[ii] = data_cleaned$foreground.log[jj]
    preProc_jointseg$replace = NA
    preProc_jointseg$replace[ii] = 1
    substitution_rate = (nrow(preProc_jointseg) - sum(is.na(preProc_jointseg$replace))) / nrow(preProc_jointseg)
    print(paste0('DryClean substitution rate: ', round(substitution_rate*100, 3), '%'))
    
    preProc_jointseg$cnlr[missing_cnlr] = NA
    preProc_jointseg$seg[missing_cnlr] = NA
    
    #' for loop alternative
    # preProc_jointseg$replace = NA
    # for(i in 1:nrow(preProc_jointseg)){
    #   if(preProc_jointseg$bin[i] %in% data_cleaned$bin){
    #     preProc_jointseg$cnlr[i] = data_cleaned$foreground.log[which(data_cleaned$bin == preProc_jointseg$bin[i])]
    #     preProc_jointseg$replace[i] = 'new'
    #   } else {
    #     preProc_jointseg$cnlr[i] = preProc_jointseg$cnlr[i]
    #     preProc_jointseg$replace[i] = 'old'
    #   }
    # }
    
    
    #' replace data frame
    preProc_jointseg = preProc_jointseg[,-c(ncol(preProc_jointseg) - 1, ncol(preProc_jointseg))]
    dat$jointseg = preProc_jointseg
    
    out = facets::procSample(dat, cval = cval, min.nhet = min_nhet)
    fit = facets::emcncf(out)
    
    # Fix bad NAs
    fit$cncf = cbind(fit$cncf, cf = out$out$cf, tcn = out$out$tcn, lcn = out$out$lcn)
    fit$cncf$lcn[fit$cncf$tcn == 1] = 0
    fit$cncf$lcn.em[fit$cncf$tcn.em == 1] = 0
    
    # Generate output
    return(list(snps = out$jointseg,
                segs = fit$cncf,
                substitution_rate = substitution_rate,
                purity = as.numeric(fit$purity),
                ploidy = as.numeric(fit$ploidy),
                dipLogR = out$dipLogR,
                alBalLogR = out$alBalLogR,
                flags = out$flags,
                em_flags = fit$emflags,
                loglik = fit$loglik))
    
    
    #' run full replacement algorithm
  } else if (MODE == 'full'){
    
    cat(paste0('Mode: ', MODE, ' is selected'))
    print(paste0('DryClean substitution rate: 100%'))
    
    # Check input 
    missing_cols = setdiff(c('Chromosome', 'Position', 'NOR.DP', 'TUM.DP', 'NOR.RD', 'TUM.RD'), names(read_counts)) 
    if (length(missing_cols) > 0) {
      stop(paste0('Input missing column(s)', paste(missing_cols, collapse = ', '), '.'), call. = FALSE)
    }
    
    set.seed(seed = seed)
    genome = 'hg19'
    
    # Run FACETS algorithm
    dat = FacetsDC::preProcSample(rcmat = read_counts, 
                                  data_cleaned = read_cleaned,
                                  ndepth = ndepth,
                                  ndepthmax = 1000,
                                  het.thresh = 0.25, 
                                  snp.nbhd = snp_nbhd, 
                                  cval = 25,
                                  gbuild = genome, 
                                  hetscale = TRUE, 
                                  unmatched = FALSE)
    
    out = facets::procSample(dat, 
                             cval = cval, 
                             min.nhet = min_nhet)
    fit = facets::emcncf(out)
    
    # Fix bad NAs
    fit$cncf = cbind(fit$cncf, cf = out$out$cf, tcn = out$out$tcn, lcn = out$out$lcn)
    fit$cncf$lcn[fit$cncf$tcn == 1] = 0
    fit$cncf$lcn.em[fit$cncf$tcn.em == 1] = 0
    
    # Generate output
    return(list(snps = out$jointseg,
                segs = fit$cncf,
                substitution_rate = '100%',
                purity = as.numeric(fit$purity),
                ploidy = as.numeric(fit$ploidy),
                dipLogR = out$dipLogR,
                alBalLogR = out$alBalLogR,
                flags = out$flags,
                em_flags = fit$emflags,
                loglik = fit$loglik))
    
  }
}
  
  