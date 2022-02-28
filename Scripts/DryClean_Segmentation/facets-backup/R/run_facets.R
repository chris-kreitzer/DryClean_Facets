## Facets compilation file:
## This file conditionally replaces CnLR estimates from Facets with those provided by DryClean
## Conditionally, because not all positions are available from DryClean

run_facets_cleaned = function(read_counts,
                              read_cleaned,
                              MODE,
                              cval = 150,
                              dipLogR = NULL,
                              ndepth = 35,
                              snp_nbhd = 250,
                              min_nhet = 15,
                              seed = 100) {
  
  set.seed(seed = seed)
  genome = 'hg19'
  
  if(missing(MODE)){
    stop("Please provide a value for MODE: ", call. = T)
  }
  
  if(!missing(MODE)){
    print(paste0('Mode: ', MODE, ' is selected'))
  
  #' Load data: count_matrix + DryClean output
  input = .readData(filename_counts = read_counts, 
                    filename_dryclean = read_cleaned)
  
  read_counts = input$rcmat
  read_cleaned = input$counts_cleaned
  print(head(read_cleaned))
  
  #' Check input 
  missing_cols = setdiff(c('Chromosome', 'Position', 'NOR.DP', 'TUM.DP', 'NOR.RD', 'TUM.RD'), names(read_counts))
  if (length(missing_cols) > 0){
    stop(paste0('Input missing column(s)', paste(missing_cols, collapse = ', '), '.'), call. = FALSE)
  }
  
  #' Run FACETS algorithm
  dat = preProcSample_DC(rcmat = read_counts,
                                data_cleaned = read_cleaned,
                                MODE = MODE,
                                ndepth = ndepth, 
                                ndepthmax = 1000,
                                het.thresh = 0.25, 
                                snp.nbhd = snp_nbhd, 
                                cval = 25,
                                gbuild = genome, 
                                hetscale = TRUE, 
                                unmatched = FALSE)
  
  out = facets::procSample(x = dat,
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
                substitution_rate = dat$substitution_rate,
                purity = as.numeric(fit$purity),
                ploidy = as.numeric(fit$ploidy),
                dipLogR = out$dipLogR,
                alBalLogR = out$alBalLogR,
                flags = out$flags,
                em_flags = fit$emflags,
                loglik = fit$loglik))
    
    
    #' run full replacement algorithm
  } else {
    
    print(paste0('Mode: ', MODE, ' is selected'))
    print(paste0('DryClean substitution rate: 100%'))
    
    # Check input 
    missing_cols = setdiff(c('Chromosome', 'Position', 'NOR.DP', 'TUM.DP', 'NOR.RD', 'TUM.RD'), names(read_counts)) 
    if (length(missing_cols) > 0) {
      stop(paste0('Input missing column(s)', paste(missing_cols, collapse = ', '), '.'), call. = FALSE)
    }
    
    # Run FACETS algorithm
    dat = FacetsDC::preProcSample(rcmat = read_counts,
                                  data_cleaned = read_cleaned,
                                  MODE = MODE,
                                  ndepth = ndepth, 
                                  ndepthmax = 1000,
                                  het.thresh = 0.25, 
                                  snp.nbhd = snp_nbhd, 
                                  cval = 25,
                                  gbuild = genome, 
                                  hetscale = TRUE, 
                                  unmatched = FALSE)
    
    out = facets::procSample(x = dat,
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
                substitution_rate = dat$substitution_rate,
                purity = as.numeric(fit$purity),
                ploidy = as.numeric(fit$ploidy),
                dipLogR = out$dipLogR,
                alBalLogR = out$alBalLogR,
                flags = out$flags,
                em_flags = fit$emflags,
                loglik = fit$loglik))
    
  }
}
  
  