##-----------------------------------------------------------------------------
## Hypersegmentation and average segment length
## ----------------------------------------------------------------------------
## 
## start: 04/23/2022
## chris-kreitzer
## 


library(DNAcopy)
clean()
gc()


MasterFile = readRDS('~/Documents/MSKCC/07_FacetsReview/DryClean/TUMOR_BRCA/tumor_table.rds')

countfiles = list.files('~/Documents/MSKCC/07_FacetsReview/DryClean/Tumor_Countfiles/', full.names = T)
cleanedfiles = list.files('~/Documents/MSKCC/07_FacetsReview/DryClean/TUMOR_DECOMPOSED/', full.names = T)


hypersegmentation = data.frame()
for(i in 1:nrow(MasterFile)){
  try({
    print(MasterFile$original[i])
    count = facetsSuite::read_snp_matrix(input_file = grep(pattern = MasterFile$original[i], x = countfiles, value = T))
    count = count[, c(1,2,3,5,4,6)]
    hyper_sub_facets = data.frame()
    for(iter in c(25, 50, 75, 100, 150, 200)){
      facets_out = facetsSuite::run_facets(read_counts = count, cval = iter, genome = 'hg19', seed = 100)
      segs = nrow(facets_out$segs)
      waterfall = dlrs(x = facets_out$snps$cnlr)
      segment_length = median(facets_out$segs$end - facets_out$segs$start)
      
      hyper_sub = data.frame(sample = MasterFile$original[i],
                             cval = iter,
                             segs = segs,
                             waterfall = waterfall,
                             median_seg_length = segment_length,
                             algorithm = 'facets')
      
      hyper_sub_facets = rbind(hyper_sub_facets, hyper_sub)
    }
    
    
    #' cleaned
    segs_clean_df_all = data.frame()
    for(cval in c(25, 50, 75, 100, 150, 200)){
      cleaned = FacetsDC::run_facets_cleaned(read_counts = grep(pattern = MasterFile$original[i], x = countfiles, value = T),
                                             read_cleaned = grep(pattern = paste0(MasterFile$sample[i], '.rds.*'), x = cleanedfiles, value = T),
                                             MODE = 'union', 
                                             cval = cval, 
                                             seed = 100)
      segs_clean = nrow(cleaned$segs)
      waterfall = dlrs(x = cleaned$snps$cnlr)
      segment_length = median(cleaned$segs$end - cleaned$segs$start)
      
      segs_clean_df = data.frame(sample = MasterFile$original[i],
                                 cval = cval,
                                 segs = segs_clean,
                                 waterfall = waterfall,
                                 median_seg_length = segment_length,
                                 algorithm = 'dryclean')
      
      segs_clean_df_all = rbind(segs_clean_df_all, segs_clean_df)
    }
    
    rm(segs_clean, waterfall, segment_length, segs_clean_df, segs, facets_out, hyper_sub)
    hypersegmentation = rbind(hypersegmentation, hyper_sub_facets, segs_clean_df_all)
  })
  
}


View(hypersegmentation)


#' Downstream
all_out = data.frame()
for(algo in unique(hypersegmentation$algorithm)){
  data = hypersegmentation[which(hypersegmentation$algorithm == algo), ]
  algo_sub = data.frame()
  for(cval in unique(data$cval)){
    median_segs = median(data$segs[which(data$cval == cval)])
    segs_CI = CI_z(data$segs[which(data$cval == cval)])
    out = data.frame(algorithm = algo,
                     cval = cval,
                     median_segs = median_segs,
                     CI = segs_CI)
    algo_sub = rbind(algo_sub, out)
    
  }
  all_out = rbind(all_out, algo_sub)
}








## CBS
facets_CBS = CNA(facets_original$snps$cnlr, facets_original$snps$chrom, facets_original$snps$maploc, data.type = 'logratio', sampleid = 'FacetsOriginal')
dryclean_CBS = CNA(cleaned$snps$cnlr, cleaned$snps$chrom, cleaned$snps$maploc, data.type = 'logratio', sampleid = 'cleaned')

smooth.facets = smooth.CNA(facets_CBS)
smooth.dryclean = smooth.CNA(dryclean_CBS)

seg_facets = segment(x = smooth.facets)
seg_dryclean = segment(x = smooth.dryclean)

par(mfrow = c(2,1))
plotSample(seg_facets, cex = 0.35, altcol = T, lwd = 2.5, pch = 20, col = c('#5c9e75', '#272425'))
plotSample(seg_dryclean, cex = 0.35, altcol = T, lwd = 2.5, pch = 20, col = c('#5c9e75', '#272425'))
dev.off()

seg_dryclean


seg_facets

plotSample(seg_facets, xlab = '', xaxt = 'n', lwd = 2, ylab = 'CBS_FacetsOriginal', main = '', las = 2)
box(lwd = 2)
dim(cleaned$segs)
dim(facets_original$segs)

a = facetsSuite::cnlr_plot(facets_original)
a1 = facetsSuite::cnlr_plot(cleaned)

a/a1

b = facetsSuite::valor_plot(facets_original)
b1 = facetsSuite::valor_plot(cleaned)

b/b1

c = facetsSuite::icn_plot(facets_original)
c1 = facetsSuite::icn_plot(cleaned)

c/c1

facets_fit_qc(facets_output = facets_original)
facets_fit_qc(facets_output = cleaned)