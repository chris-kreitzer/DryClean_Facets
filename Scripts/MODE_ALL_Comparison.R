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
for(i in 1:length(Facets_files)){
  name = basename(Facets_files[i])
  name = substr(name, start = 17, stop = 47)
  print(name)
  data_in = facetsSuite::read_snp_matrix(input_file = Facets_files[i])
  data_processed = facetsSuite::run_facets(read_counts = data_in, cval = 150, ndepth = 35, snp_nbhd = 250, seed = 100)
  purity = data_processed$purity
  ploidy = data_processed$ploidy
  fit = data_processed$segs
  fit$name = name
  QC = facets_fit_qc(data_processed)$facets_qc
  
  summary = data.frame(sample = name, 
                       purity = purity,
                       ploidy = ploidy,
                       n.segs = nrow(fit),
                       QC = QC)
  
  Facets_original_segments = rbind(Facets_original_segments, fit)
  Facets_original_summary = rbind(Facets_original_summary, summary)
  
  rm(data_in, data_processed, purity, ploidy, QC, summary, fit)
  
}

Facets_original_summary$sample =  substr(Facets_original_summary$sample, start = 17, stop = 47)


#' Facets MODE1 IMPLEMENATION (Partial replacement):
## run updated Facets on the samples
Dryclean_files = list.files('Tumor_cleaned/', full.names = T)

DryClean_segments = data.frame()
DryClean_summary = data.frame()
DryClean_QC = data.frame()
for(i in 1:length(Facets_files)){
  name = basename(Facets_files[i])
  name = substr(name, start = 17, stop = 47)
  print(name)
  name.dry = grep(pattern = name, x = Dryclean_files, value = T)
  data_in = FacetsDC::readSnpMatrix(filename_counts = Facets_files[i],
                                    filename_dryclean = name.dry)
  data_out = FacetsDC::run_facets_cleaned(read_counts = data_in$rcmat,
                                read_cleaned = data_in$counts_cleaned, 
                                MODE = 'partial', 
                                cval = 150, 
                                ndepth = 35, 
                                snp_nbhd = 250, 
                                seed = 100)
  purity = data_out$purity
  ploidy = data_out$ploidy
  fit = data_out$segs
  fit$name = name
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
  rm(data_in, data_out, purity, ploidy, fit, QC, QC_all)
}


#' Facets MODE2 IMPLEMENATION (Full replacement, information loss):
## run updated Facets on the samples
Facets_files
Dryclean_files = list.files('Tumor_cleaned/', full.names = T)

DryClean_segments_full = data.frame()
DryClean_summary_full = data.frame()
DryClean_QC_full = data.frame()
for(i in 1:length(Facets_files)){
  name = basename(Facets_files[i])
  name = substr(name, start = 17, stop = 47)
  print(name)
  name.dry = grep(pattern = name, x = Dryclean_files, value = T)
  data_in = FacetsDC::readSnpMatrix(filename_counts = Facets_files[i],
                                    filename_dryclean = name.dry)
  data_out = run_facets_cleaned(read_counts = data_in$rcmat,
                                read_cleaned = data_in$counts_cleaned, 
                                MODE = 'full', 
                                cval = 150, 
                                ndepth = 35, 
                                snp_nbhd = 250, 
                                seed = 100)
  purity = data_out$purity
  ploidy = data_out$ploidy
  fit = data_out$segs
  fit$name = name
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
  rm(data_in, data_out, purity, ploidy, fit, QC, QC_all)
}


#' summary stats
df_list = cbind(Facets_original_summary, DryClean_summary, DryClean_summary_full)

