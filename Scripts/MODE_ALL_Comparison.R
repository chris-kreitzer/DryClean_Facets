#' Comparing the two modes of Implemenation:

rm(list = ls())
.rs.restartR()
setwd('~/Documents/GitHub/DryClean_Facets/')
source('Scripts/FacetsQC.R')
library(FacetsDC)
library(facets)

## run Facets on all files
Facets_files = list.files(path = 'Tumor_countsFile/', full.names = T)
set.seed(100)
Facets_original_segments = data.frame()
Facets_original_summary = data.frame()
Facets_original_QC = data.frame()
Facets_original_snps = data.frame()
for(i in 64:length(Facets_files)){
  name = basename(Facets_files[i])
  name = substr(name, start = 17, stop = 47)
  print(name)
  data_in = facetsSuite::read_snp_matrix(input_file = Facets_files[i])
  data_pre = facets::preProcSample(rcmat = data_in, 
                                   het.thresh = 0.25, 
                                   snp.nbhd = 250, 
                                   ndepthmax = 1000, 
                                   cval = 25, 
                                   ndepth = 35)
  data_post = facets::procSample(data_pre, cval = 150)
  data_fit = facets::emcncf(data_post)
  fit = data_fit$cncf
  fit$name = name
  purity = data_fit$purity
  ploidy = data_fit$ploidy
  snps = data_post$jointseg
  snps$name = name
  
  # QC = facets_fit_qc(data_processed)$facets_qc
  # QC_all = facets_fit_qc(facets_output = data_processed)
  # QC_all$name = name
  
  summary = data.frame(sample = name, 
                       purity = purity,
                       ploidy = ploidy,
                       n.segs = nrow(fit))
  
  Facets_original_segments = rbind(Facets_original_segments, fit)
  Facets_original_summary = rbind(Facets_original_summary, summary)
  #Facets_original_QC = rbind(Facets_original_QC, QC_all)
  Facets_original_snps = rbind(Facets_original_snps, snps)
  
  rm(data_in, data_processed, purity, ploidy, QC, summary, fit, QC_all, snps, data_in)
  
}

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
FD = merge(Facets_original_summary, DryClean_summary, by = 'sample')
FDE = merge(FD, DryClean_summary_full, by = 'sample')
colnames(FDE) = c('sample', 'Facets_purity', 'Facets_ploidy', 'facets_n.segs', 
                  'facets_QC', 'Dryclean_purity', 'Dryclean_ploidy', 
                  'DryClean_segs', 'Dryclean_subst.', 'Dryclean_QC', 
                  'DCF_purity', 'DCF_ploidy', 'DCF_segs', 'DCF_subst', 'DCF_QC')


#' A potential information loss:
ggplot(FDE, aes(x = facets_n.segs, y = DCF_segs)) + 
  geom_jitter() +
  scale_y_continuous(limits = c(20,120)) +
  scale_x_continuous(limits = c(20,120))+
  geom_abline(slope = 1, intercept = 0) +
  theme_minimal() +
  theme(aspect.ratio = 1,
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 12, color = 'black')) +
  labs(x = '# segments Facets', y = '# segments DryClean')


comparision_df_F = Facets_original_summary
comparision_df_F$type = rep('Facets', nrow(comparision_df_F))
comparision_df_DF = DryClean_summary_full
comparision_df_DF$type = rep('DryClean', nrow(comparision_df_DF))
comparision_df_DF$substitution_rate = NULL
comparision = rbind(comparision_df_F, comparision_df_DF)
View(comparision)

ggplot(comparision, aes(x = type, y = purity)) + 
  geom_boxplot(width = 0.6) +
  geom_jitter(shape = 16, position = position_jitter(0.1)) +
  scale_y_continuous(limits = c(0,1)) +
  theme_minimal() +
  theme(aspect.ratio = 2,
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 12, color = 'black'),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(x = '', y = 'Purity', title = 'Full replacement')






#' Replacement of estimates completely random - follow certain pattern or uniform?
#' sample with lowest substitution rate: P-0038760-T01-IM6_P-0038760-N01 (56 %)

data_in = FacetsDC::readSnpMatrix(filename_counts = 'Tumor_countsFile/countsMerged____P-0038760-T01-IM6_P-0038760-N01-IM6.dat.gz',
                                  filename_dryclean = 'Tumor_cleaned/P-0038760-T01-IM6_P-0038760-N01-IM6_drycleaned.rds')

set.seed(100)
fo = facets::preProcSample(rcmat = data_in$rcmat,
                           het.thresh = 0.25, 
                           snp.nbhd = 250, 
                           ndepthmax = 1000, 
                           cval = 25, 
                           ndepth = 35)
fo = facets::procSample(fo, cval = 150)
fo_joint = fo$jointseg
fo_joint$bin = paste(fo_joint$chrom, fo_joint$maploc, sep = ';')
head(fo_joint)


dm = FacetsDC::run_facets_cleaned(read_counts = data_in$rcmat, 
                                  read_cleaned = data_in$counts_cleaned,
                                  MODE = 'partial', 
                                  cval = 150, 
                                  snp_nbhd = 250, 
                                  seed = 100)
dm_joint = dm$snps
dm_joint$bin = paste(dm_joint$chrom, dm_joint$maploc, sep = ';')

head(dm_joint)

setdiff(fo_joint$bin, dm_joint$bin)
setdiff(dm_joint$bin, fo_joint$bin)

length(intersect(dm_joint$bin, fo_joint$bin))
head(dm_joint)
head(fo_joint)



View(dm_joint)
View(fo_joint)



















