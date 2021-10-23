## Modification of Facets algorithm; 
## substitute the raw (counts) data with the cleaned input from DryClean;
## retain (Facets-) Positions where no DryClean estimate is available == MODE 1
## 
## 10/20/21
## chris-kreitzer
## 


#' downstream
setwd('~/Documents/GitHub/DryClean_Facets/')

## run Facets on all files
Facets_files = list.files(path = 'Tumor_countsFile/', full.names = T)

Facets_original_segments = data.frame()
Facets_original_summary = data.frame()
for(i in 1:length(Facets_files)){
  name = basename(Facets_files[i])
  print(name)
  data_in = facetsSuite::read_snp_matrix(input_file = Facets_files[i])
  data_processed = facetsSuite::run_facets(read_counts = data_in)
  purity = data_processed$purity
  ploidy = data_processed$ploidy
  fit = data_processed$segs
  fit$name = name
  
  summary = data.frame(sample = name, 
                       purity = purity,
                       ploidy = ploidy,
                       n.segs = nrow(fit))
  
  Facets_original_segments = rbind(Facets_original_segments, fit)
  Facets_original_summary = rbind(Facets_original_summary, summary)
  
  rm(data_in, data_processed, purity, ploidy)
  
}


## run updated Facets on the samples
Facets_files
Dryclean_files = list.files('Tumor_cleaned/', full.names = T)

DryClean_segments = data.frame()
DryClean_summary = data.frame()
for(i in 1:length(Facets_files)){
  name = basename(Facets_files[i])
  name = substr(name, start = 17, stop = 47)
  print(name)
  name.dry = grep(pattern = name, x = Dryclean_files, value = T)
  data_in = FacetsDC::readSnpMatrix(filename_counts = Facets_files[i],
                                    filename_dryclean = name.dry)
  data_out = FacetsDC::run_facets_cleaned(read_counts = data_in$rcmat,
                                          read_cleaned = data_in$counts_cleaned)
  purity = data_out$purity
  ploidy = data_out$ploidy
  fit = data_out$segs
  fit$name = name
  
  summary_dry = data.frame(sample = name,
                           purity = purity,
                           ploidy = ploidy,
                           n.segs = nrow(fit))
  
  DryClean_segments = rbind(DryClean_segments, fit)
  DryClean_summary = rbind(DryClean_summary, summary_dry)
  rm(data_in, data_out, purity, ploidy, fit)
}




#' merge the data and compare:
Facets_original_summary$sample = substr(Facets_original_summary$sample, start = 17, stop = 47)
MODE1 = merge(Facets_original_summary, DryClean_summary, by = 'sample')

#' investiage the output
library(ggplot2)
n_segs = ggplot(MODE1, aes(x = n.segs.x, y = n.segs.y)) + 
  geom_point() +
  geom_abline(slope = 1, linetype = "dashed", color = 'grey35', size = 0.2) +
  scale_y_continuous(limits = c(0,120), expand = c(0,0)) +
  scale_x_continuous(limits = c(0,120), expand = c(0,0)) +
  theme_classic() +
  theme(aspect.ratio = 1,
        panel.border = element_rect(fill = NA, size = 0.8)) +
  labs(title = 'MODE1: all sampled Facets-loci;\nloci substituted where applicable', x = '# segments Facets', y = '# segments modified Facets')
n_segs

ggsave_golden(n_segs, filename = 'Figures/MODE1_n.segs.pdf', width = 8)


#' split into low and high purity samples:
BRCA = readRDS('Data_out/BRCA_workingCohort_MSK.rds')
BRCA_high = BRCA$tumor_sample[which(BRCA$purity >= 0.65)]

high_purity = data.frame()
for(i in 1:length(BRCA_high)){
  out = MODE1[grepl(pattern = BRCA_high[i], x = MODE1$sample), ]
  high_purity = rbind(high_purity, out)
}

n_segs_high = ggplot(high_purity, aes(x = n.segs.x, y = n.segs.y)) + 
  geom_point() +
  geom_abline(slope = 1, linetype = "dashed", color = 'grey35', size = 0.2) +
  scale_y_continuous(limits = c(0,120), expand = c(0,0)) +
  scale_x_continuous(limits = c(0,120), expand = c(0,0)) +
  theme_classic() +
  theme(aspect.ratio = 1,
        panel.border = element_rect(fill = NA, size = 0.8)) +
  labs(title = 'High purity samples [0.65 < x < 1]', x = '# segments Facets', y = '# segments modified Facets')

n_segs_high


#' low purity samples:
BRCA_low = BRCA$tumor_sample[which(BRCA$purity < 0.65)]

low_purity = data.frame()
for(i in 1:length(BRCA_low)){
  out = MODE1[grepl(pattern = BRCA_low[i], x = MODE1$sample), ]
  low_purity = rbind(low_purity, out)
}

n_segs_low = ggplot(low_purity, aes(x = n.segs.x, y = n.segs.y)) + 
  geom_point() +
  geom_abline(slope = 1, linetype = "dashed", color = 'grey35', size = 0.2) +
  scale_y_continuous(limits = c(0,120), expand = c(0,0)) +
  scale_x_continuous(limits = c(0,120), expand = c(0,0)) +
  theme_classic() +
  theme(aspect.ratio = 1,
        panel.border = element_rect(fill = NA, size = 0.8)) +
  labs(title = 'Low purity samples [0.0 < x < 0.65]', x = '# segments Facets', y = '# segments modified Facets')

n_segs_low

library(patchwork)

MODE1_out = n_segs + (n_segs_high / n_segs_low)
ggsave_golden(filename = 'Figures/MODE1_comparision.pdf', plot = MODE1_out, width = 12)




#' breakpoints:
Facets_original_segments$length = (Facets_original_segments$end - Facets_original_segments$start) / 1000000
Facets_original_segments$rel_marker = Facets_original_segments$num.mark / Facets_original_segments$length

ggplot(Facets_original_segments[which(Facets_original_segments$rel_marker < 100), ], aes(x = rel_marker)) +
  geom_density(aes(y = ..density..), colour = "black", fill = "white", size = 1.2) +
  scale_y_continuous(expand = c(0.05, 0)) + 
  scale_x_continuous(expand = c(0.01, 0.01)) +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        panel.border = element_rect(fill = NA)) +
  labs(title = 'Facets: Markers per Mb', x = '# markers', y = 'density')


DryClean_segments$length = (DryClean_segments$end - DryClean_segments$start) / 1000000
DryClean_segments$rel_marker = DryClean_segments$num.mark / DryClean_segments$length

ggplot(DryClean_segments[which(DryClean_segments$rel_marker < 100), ], aes(x = rel_marker)) +
  geom_density(aes(y = ..density..), colour = "black", fill = "white", size = 1.2) +
  scale_y_continuous(expand = c(0.05, 0)) + 
  scale_x_continuous(expand = c(0.01, 0.01)) +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        panel.border = element_rect(fill = NA)) +
  labs(title = 'DryClean: Markers per Mb', x = '# markers', y = 'density')

median(DryClean_segments$rel_marker)
median(Facets_original_segments$rel_marker)
