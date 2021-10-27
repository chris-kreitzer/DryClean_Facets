#' Comparing the two modes of Implemenation:

rm(list = ls())
.rs.restartR()
setwd('~/Documents/GitHub/DryClean_Facets/')
source('Scripts/FacetsQC.R')
library(FacetsDC)
library(facets)

## run Facets on all files
Facets_files = list.files(path = 'Tumor_countsFile/', full.names = T)

Facets_original_segments = data.frame()
Facets_original_summary = data.frame()
Facets_original_QC = data.frame()
Facets_original_snps = data.frame()
for(i in 1:length(Facets_files)){
  name = basename(Facets_files[i])
  name = substr(name, start = 17, stop = 47)
  print(name)
  data_in = facetsSuite::read_snp_matrix(input_file = Facets_files[i])
  data_fit = facetsSuite::run_facets(read_counts = data_in, 
                                     cval = 150, 
                                     ndepth = 35,
                                     dipLogR = NULL, 
                                     snp_nbhd = 250, 
                                     min_nhet = 15, 
                                     genome = 'hg19', 
                                     seed = 100) 
                                   
  fit = data_fit$segs
  fit$name = name
  purity = data_fit$purity
  ploidy = data_fit$ploidy
  snps = data_fit$snps
  snps$name = name
  
  QC = facets_fit_qc(data_fit)$facets_qc
  QC_all = facets_fit_qc(facets_output = data_fit)
  QC_all$name = name
  
  summary = data.frame(sample = name, 
                       purity = purity,
                       ploidy = ploidy,
                       n.segs = nrow(fit),
                       QC = QC)
  
  Facets_original_segments = rbind(Facets_original_segments, fit)
  Facets_original_summary = rbind(Facets_original_summary, summary)
  Facets_original_QC = rbind(Facets_original_QC, QC_all)
  Facets_original_snps = rbind(Facets_original_snps, snps)
  
  rm(data_in, data_processed, purity, ploidy, QC, summary, fit, QC_all, snps, data_in)
  
}

#' Facets MODE1 IMPLEMENATION (Partial replacement):
## run updated Facets on the samples
Dryclean_files = list.files('Tumor_cleaned/', full.names = T)

DryClean_segments = data.frame()
DryClean_summary = data.frame()
DryClean_QC = data.frame()
DryClean_snps = data.frame()
for(i in 1:length(Facets_files)){
  name = basename(Facets_files[i])
  name = substr(name, start = 17, stop = 47)
  print(name)
  name.dry = grep(pattern = name, x = Dryclean_files, value = T)
  data_in = FacetsDC::readSnpMatrix(filename_counts = Facets_files[i],
                                    filename_dryclean = name.dry)
  data_out = FacetsDC::run_facets_cleaned(read_counts = data_in$rcmat, 
                                          read_cleaned = data_in$counts_cleaned, 
                                          dipLogR = NULL, 
                                          min_nhet = 15,
                                          MODE = 'partial', 
                                          cval = 150, 
                                          seed = 100, 
                                          snp_nbhd = 250)
  purity = data_out$purity
  ploidy = data_out$ploidy
  fit = data_out$segs
  fit$name = name
  snps = data_out$snps
  snps$name = name
  substitution_rate = data_out$substitution_rate
  QC = facets_fit_qc(facets_output = data_out)$facets_qc
  QC_all = facets_fit_qc(facets_output = data_out)
  
  summary_dry = data.frame(sample = name,
                           purity = purity,
                           ploidy = ploidy,
                           n.segs = nrow(fit),
                           substitution_rate = substitution_rate,
                           QC = QC)
  
  DryClean_segments = rbind(DryClean_segments, fit)
  DryClean_summary = rbind(DryClean_summary, summary_dry)
  DryClean_QC = rbind(DryClean_QC, QC_all)
  DryClean_snps = rbind(DryClean_snps, snps)
  rm(data_in, data_out, purity, ploidy, fit, QC, QC_all)
}


#' Facets MODE2 IMPLEMENATION (Full replacement, information loss):
## run updated Facets on the samples
Facets_files
Dryclean_files = list.files('Tumor_cleaned/', full.names = T)

DryClean_segments_full = data.frame()
DryClean_summary_full = data.frame()
DryClean_QC_full = data.frame()
DryClean_snps_full = data.frame()
for(i in 1:length(Facets_files)){
  name = basename(Facets_files[i])
  name = substr(name, start = 17, stop = 47)
  print(name)
  name.dry = grep(pattern = name, x = Dryclean_files, value = T)
  data_in = FacetsDC::readSnpMatrix(filename_counts = Facets_files[i],
                                    filename_dryclean = name.dry)
  
  data_out = FacetsDC::run_facets_cleaned(read_counts = data_in$rcmat, 
                                          read_cleaned = data_in$counts_cleaned, 
                                          dipLogR = NULL, 
                                          min_nhet = 15,
                                          MODE = 'full', 
                                          cval = 150, 
                                          seed = 100, 
                                          snp_nbhd = 250)
  purity = data_out$purity
  ploidy = data_out$ploidy
  fit = data_out$segs
  fit$name = name
  snps = data_out$snps
  snps$name = name
  substitution_rate = data_out$substitution_rate
  QC = facets_fit_qc(facets_output = data_out)$facets_qc
  QC_all = facets_fit_qc(facets_output = data_out)
  
  summary_dry = data.frame(sample = name,
                           purity = purity,
                           ploidy = ploidy,
                           n.segs = nrow(fit),
                           substitution_rate = substitution_rate,
                           QC = QC)
  
  DryClean_segments_full = rbind(DryClean_segments_full, fit)
  DryClean_summary_full = rbind(DryClean_summary_full, summary_dry)
  DryClean_QC_full = rbind(DryClean_QC_full, QC_all)
  DryClean_snps_full = rbind(DryClean_snps_full, snps)
  rm(data_in, data_out, purity, ploidy, fit, QC, QC_all)
}


#' summary stats
FD = merge(Facets_original_summary, DryClean_summary, by = 'sample')
FDE = merge(FD, DryClean_summary_full, by = 'sample')
colnames(FDE) = c('sample', 'Facets_purity', 'Facets_ploidy', 'facets_n.segs', 
                  'Dryclean_purity', 'Dryclean_ploidy', 
                  'DryClean_segs', 'Dryclean_subst.', 'Dryclean_QC', 
                  'DCF_purity', 'DCF_ploidy', 'DCF_segs', 'DCF_subst', 'DCF_QC')


#' A potential information loss:
partial = ggplot(FDE, aes(x = facets_n.segs, y = DryClean_segs)) + 
  geom_jitter() +
  scale_y_continuous(limits = c(20, 100)) +
  scale_x_continuous(limits = c(20, 100))+
  geom_abline(slope = 1, intercept = 0) +
  theme_minimal() +
  theme(aspect.ratio = 1,
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 12, color = 'black')) +
  labs(x = '# segments Facets', y = '# segments DryClean [partial]')

full = ggplot(FDE, aes(x = facets_n.segs, y = DCF_segs)) + 
  geom_jitter() +
  scale_y_continuous(limits = c(20, 100)) +
  scale_x_continuous(limits = c(20, 100))+
  geom_abline(slope = 1, intercept = 0) +
  theme_minimal() +
  theme(aspect.ratio = 1,
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 12, color = 'black')) +
  labs(x = '# segments Facets', y = '# segments DryClean [partial]')




#





