## parameters 
#' ndepth = 35
#' het.thresh = 0.25
#' ndepthmax = 1000
# are neglected, because those parameters are already considered in the DryClean pipeline

procSnps = function(rcmat, nX = 23) {
    # keep only chromsomes 1-22 & X for humans
    # modify the rcmat input data frame
    rcmat = rcmat[, c('seqnames', 'start', 'foreground.log')]
    colnames(rcmat) = c('Chromsome', 'Position', 'foreground.log')
    
    chromlevels = c(1:(nX-1), "X")
    chr.keep = rcmat$Chromosome %in% chromlevels
    
    # output data frame
    out = list()
    out$chrom = rcmat$Chromosome
    out$maploc = rcmat$Position
    out$cnlr = rcmat$foreground.log
    out$valor = 0
    out$het = 0
    #out$rCountT <- rcmat$TUM.DP
    #out$rCountN <- rcmat$NOR.DP
    
    #' set the variant allele frequency to zero, as there is no allele-specific stream yet
    # out$vafT = 0
    # out$vafN = 0
    
    # make chromosome ordered and numeric
    out$chrom = as.numeric(ordered(out$chrom, levels = chromlevels))
    
    # # call a snp heterozygous if min(vafN, 1-mafN) > het.thresh
    # if (unmatched) {
    #     if (het.thresh == 0.25) het.thresh <- 0.1
    #     out$het <- 1*(pmin(out$vafT, 1-out$vafT) > het.thresh & out$rCountT >= 50)
    # } else {
    #     out$het <- 1*(pmin(out$vafN, 1-out$vafN) > het.thresh)
    # }
    # scan maploc for snps that are close to one another (within snp.nbhd bases)
    # heep all the hets (should change if too close) and only one from a nbhd
    # out$keep = scanSnp(out$maploc, out$het, snp.nbhd)
    as.data.frame(out)
    message(print(head(out)))
}

# scanSnp <- function(maploc, het, nbhd) {
#     n <- length(maploc)
#     zzz <- .Fortran("scansnp",
#                     as.integer(n),
#                     as.double(maploc),
#                     as.double(het),
#                     keep=double(n),
#                     as.double(nbhd))
#     zzz$keep
# }

# obtain logR and logOR from read counts and GC-correct logR
counts2logROR = function(mat, gbuild, ugcpct = NULL, f = 0.2) {
    out = mat
    # gc percentage
    out$gcpct <- rep(NA_real_, nrow(out))
    # get GC percentages from pctGCdata package
    # loop thru chromosomes
    nchr <- max(mat$chrom) # IMPACT doesn't have X so only 22
    for (i in 1:nchr) {
        ii <- which(out$chrom==i)
        # allow for chromosomes with no SNPs i.e. not targeted
        if (length(ii) > 0) {
            if (gbuild == "udef") {
                out$gcpct[ii] <- getGCpct(i, out$maploc[ii], gbuild, ugcpct)
            } else {
                out$gcpct[ii] <- getGCpct(i, out$maploc[ii], gbuild)
            }
        }
    }
    ##### log-ratio with gc correction and maf log-odds ratio steps
    chrom = out$chrom
    maploc = out$maploc
    rCountN = 0
    rCountT = 0
    vafT = 0
    vafN = 0
    het = 0
    gcpct = 0
    
    # compute gc bias
    # ncount <- tapply(rCountN, gcpct, sum)
    # tcount <- tapply(rCountT, gcpct, sum)
    # pctgc <- as.numeric(names(ncount))
    # tscl <- sum(ncount)/sum(tcount)
    # gcb <- lowess(pctgc, log2(tcount*tscl)-log2(ncount), f=f)
    # jj <- match(gcpct, gcb$x)
    # gcbias <- gcb$y[jj]
    
    
    # compute cn log-ratio (gc corrected) and baf log odds-ratio
    #cnlr <- log2(1+rCountT*tscl) - log2(1+rCountN) - gcbias
    cnlr = out$foreground.log
    
    
    
    # minor allele log-odds ratio and weights
    #lorvar = valor <- rep(NA_real_, length(maploc))
    lorvar = 0
    # if (unmatched) {
    #     # read count matrix for odds ratio etc
    #     rcmat <- round(cbind(vafT[het==1]*rCountT[het==1], (1-vafT[het==1])*rCountT[het==1]))
    #     # folded log of Tukey (with 1/6 correction)
    #     valor[het==1] <- log(rcmat[,1]+1/6) - log(rcmat[,2]+1/6)
    #     # variance - approximation using delta method
    #     lorvar[het==1] <- 1/(rcmat[,1]+1/6) + 1/(rcmat[,2]+1/6)
    # } else {
    #     # read count matrix for odds ratio etc
    #     rcmat <- round(cbind(vafT[het==1]*rCountT[het==1], (1-vafT[het==1])*rCountT[het==1], vafN[het==1]*rCountN[het==1], (1-vafN[het==1])*rCountN[het==1]))
    #     # log-odds-ratio (Haldane correction)
    #     valor[het==1] <- log(rcmat[,1]+0.5) - log(rcmat[,2]+0.5) - log(rcmat[,3]+0.5) + log(rcmat[,4]+0.5)
    #     # variance of log-odds-ratio (Haldane; Gart & Zweifel Biometrika 1967)
    #     lorvar[het==1] <- (1/(rcmat[,1]+0.5) + 1/(rcmat[,2]+0.5) + 1/(rcmat[,3]+0.5) + 1/(rcmat[,4]+0.5))
    # }
    
    # put them together
    out$lorvar = 0
    out$gcbias = 0
    out$cnlr = cnlr
    out$valor = 0
    out$lorvar = 0
    out$het = 0
    out
}


