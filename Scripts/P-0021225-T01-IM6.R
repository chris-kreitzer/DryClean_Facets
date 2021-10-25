## Case Study: P-0021225-T01-IM6; primary breast cancer patient (age = 55)
## Note: Low tumor content (approximately 20% or less). Negative mutation, 
## copy number and structural variant results, as well as low (stable) MSIsensor scores, 
## should be interpreted with caution (cBIO)
## 
## potential issue: missing true positives


setwd('~/Documents/GitHub/DryClean_Facets/')
cohort = readRDS('Data_out/BRCA_workingCohort_MSK.rds')
cohort$tumor_sample[which(cohort$purity < 0.2)]
library(FacetsDC)
library(facets)


#' First example
sample = 'P-0021225-T01-IM6'


#' run the pipeline:
Paths = list.files('Tumor_countsFile/', full.names = T)
Paths_cleaned = list.files('Tumor_cleaned/', full.names = T)
select_sample = which(grepl(pattern = sample, x = Paths))
select_sample_cleaned = which(grepl(pattern = sample, x = Paths_cleaned))

dat = facetsSuite::read_snp_matrix(input_file = Paths[select_sample])
dat_fit = facetsSuite::run_facets(read_counts = dat,
                                  cval = 150, 
                                  snp_nbhd = 250,
                                  genome = 'hg19')
facetsSuite::icn_plot(dat_fit)


#' run DryClean: MODE1
dat_dryclean = FacetsDC::readSnpMatrix(filename_counts = Paths[select_sample],
                                             filename_dryclean = Paths_cleaned[select_sample_cleaned])

dat_dryclean_mode1_fit = FacetsDC::run_facets_cleaned(read_counts = dat_dryclean_mode1$rcmat,
                                                      read_cleaned = dat_dryclean_mode1$counts_cleaned, 
                                                      MODE = 'partial', cval = 150, snp_nbhd = 250)
facetsSuite::icn_plot(dat_dryclean_mode1_fit)


#' run DryClean: MODE2
dat_dryclean_mode2_fit = FacetsDC::run_facets_cleaned(read_counts = dat_dryclean_mode1$rcmat,
                                                      read_cleaned = dat_dryclean_mode1$counts_cleaned, 
                                                      MODE = 'full', cval = 150, snp_nbhd = 250)
dat_dryclean_mode2_fit$segs


#' Investiage the signals at chromosome 3:
DryClean = dat_dryclean$counts_cleaned
DryClean_chr3 = DryClean[which(DryClean$seqnames == 3), ]

foreground = ggplot(DryClean_chr3, aes(x = start, y = foreground.log)) +
  geom_point() + 
  scale_y_continuous(expand = c(0, 0), 
                     limits = c(-1.5, 1.5)) +
  scale_x_continuous(expand = c(0,0)) +
  geom_hline(yintercept = seq(-1.5, 1.5, by = 0.5),
             linetype = 'dashed', size = 0.2) +
  geom_segment(aes(x = PIK3CA.loc[1], xend = PIK3CA.loc[length(PIK3CA.loc)], y = -1.2, yend = 1.2),
               size = 0.1, color = 'red', linetype = 'dashed') +
  theme(axis.line.y = element_line(size = 0.7, color = 'black'),
        axis.ticks.y = element_line(size = 0.7, lineend = 'round', color = 'black'), 
        axis.ticks.length.y = unit(0.2, 'cm'),
        panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 10, color = 'black'),
        axis.text.x = element_blank()) +
  annotate("rect", 
           xmin = -10, 
           xmax = max(DryClean_chr3$start), 
           ymin = -1.5, 
           ymax = 1.5, 
           alpha = 0.1,
           color = 'grey') +
  labs(x = '', y = 'Log foreground')

background = ggplot(DryClean_chr3, aes(x = start, y = background.log)) +
  geom_point() + 
  scale_y_continuous(expand = c(0, 0), 
                     limits = c(-1.5, 1.5)) +
  scale_x_continuous(expand = c(0,0)) +
  geom_hline(yintercept = seq(-1.5, 1.5, by = 0.5),
             linetype = 'dashed', size = 0.2) +
  geom_segment(aes(x = PIK3CA.loc[1], xend = PIK3CA.loc[length(PIK3CA.loc)], y = -1.2, yend = 1.2),
               size = 0.1, color = 'red', linetype = 'dashed') +
  theme(axis.line.y = element_line(size = 0.7, color = 'black'),
        axis.ticks.y = element_line(size = 0.7, lineend = 'round', color = 'black'), 
        axis.ticks.length.y = unit(0.2, 'cm'),
        panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 10, color = 'black'),
        axis.text.x = element_blank()) +
  annotate("rect", 
           xmin = -10, 
           xmax = max(DryClean_chr3$start), 
           ymin = -1.5, 
           ymax = 1.5, 
           alpha = 0.1,
           color = 'grey') +
  labs(x = '', y = 'Log background')


#' Facets full mode
Facets_original = dat_fit$snps[which(dat_fit$snps$chrom == 3), ]
Facets_original$replace = NA
ii = which(Facets_original$maploc %in% DryClean_chr3$start)
Facets_original$replace[ii] = 'replaced_Dryclean'
Facets_original$replace[which(is.na(Facets_original$replace))] = 'original_cnlr'

original = ggplot(Facets_original, aes(x = maploc, y = cnlr, color = replace)) +
  geom_point() +
  scale_y_continuous(expand = c(0, 0), 
                     limits = c(-1.5, 1.5)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_color_manual(values = c('original_cnlr' = 'red',
                                'replaced_Dryclean' = 'blue')) +
  geom_hline(yintercept = seq(-1.5, 1.5, by = 0.5),
             linetype = 'dashed', size = 0.2) +
  geom_segment(aes(x = PIK3CA.loc[1], xend = PIK3CA.loc[length(PIK3CA.loc)], y = -1.2, yend = 1.2),
               size = 0.1, color = 'red', linetype = 'dashed') +
  theme(axis.line.y = element_line(size = 0.7, color = 'black'),
        axis.ticks.y = element_line(size = 0.7), 
        axis.ticks.length.y = unit(0.2, 'cm'),
        panel.background = element_blank(),
        axis.text.y = element_text(size = 10, color = 'black'),
        axis.text.x = element_blank(),
        axis.line.x = element_line(color = 'black', size = 0.7),
        axis.ticks.length.x = unit(0.2, 'cm')) +
  annotate("rect", 
           xmin = -10, 
           xmax = max(DryClean_chr3$start), 
           ymin = -1.5, 
           ymax = 1.5, 
           alpha = 0.1,
           color = 'grey') +
  labs(x = 'Chr. 3', y = 'CnLR [Facets/DryClean]')



ggsave_golden(filename = 'Figures/P-0021225-T01-IM6_PIK3CA.pdf', plot = (foreground / background / original), width = 12)


#' mark PIK3CA
start = 178866145
end = 178957881
PIK3CA.loc = DryClean_chr3$start[which(DryClean_chr3$start >= start & DryClean_chr3$start <= end)]


#' compare the PIK3CA locus from DryClean and Facets:
Facets_PIK3CA = dat_fit$snps[which(dat_fit$snps$chrom == 3 & dat_fit$snps$maploc >= start & dat_fit$snps$maploc <= end), ]
DryClean_PIK3CA = DryClean_chr3[which(DryClean_chr3$start >= start & DryClean_chr3$start <= end), ]
TN_PIK3CA = dat[which(dat$Chromosome == 3 & dat$Position >= start & dat$Position <= end), ]
TN_PIK3CA$ratio = log(TN_PIK3CA$NOR.DP / TN_PIK3CA$TUM.DP)
TN_PIK3CA$ratio[is.infinite(TN_PIK3CA$ratio)] = 0

DryClean_PIK3CA = ggplot(DryClean_PIK3CA, aes(x = start, y = foreground.log)) +
  geom_point() +
  scale_y_continuous(expand = c(0, 0), 
                     limits = c(-1.5, 1.5)) +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  geom_hline(yintercept = seq(-1.5, 1.5, by = 0.5),
             linetype = 'dashed', size = 0.2) +
  theme(axis.line.y = element_line(size = 0.7, color = 'black'),
        axis.ticks.y = element_line(size = 0.7, lineend = 'round', color = 'black'), 
        axis.ticks.length.y = unit(0.2, 'cm'),
        panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 14, color = 'black'),
        axis.title.y = element_text(size = 14, color = 'black'),
        axis.text.x = element_blank(),
        aspect.ratio = 1) +
  annotate('text',
           x = 178937025,
           y = 1.2,
           size = 8,
           label = paste0('dlrs: ', round(dlrs(x = DryClean_PIK3CA$foreground.log), 3))) +
  labs(x = '', y = 'Log foreground')

Facets_PIK3CA = ggplot(Facets_PIK3CA, aes(x = maploc, y = cnlr)) +
  geom_point() +
  scale_y_continuous(expand = c(0, 0), 
                     limits = c(-1.5, 1.5)) +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  geom_hline(yintercept = seq(-1.5, 1.5, by = 0.5),
             linetype = 'dashed', size = 0.2) +
  theme(axis.line.y = element_line(size = 0.7, color = 'black'),
        axis.ticks.y = element_line(size = 0.7, lineend = 'round', color = 'black'), 
        axis.ticks.length.y = unit(0.2, 'cm'),
        panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 14, color = 'black'),
        axis.text.x = element_blank(),
        axis.line.x = element_line(size = 0.7, color = 'black'),
        axis.title.y = element_text(size = 14, color = 'black'),
        aspect.ratio = 1) +
  annotate('text',
           x = 178937025,
           y = 1.2,
           size = 8,
           label = paste0('dlrs: ', round(dlrs(x = Facets_PIK3CA$cnlr), 3))) +
  labs(x = 'PIK3CA locus', y = 'CnLR')

#' raw TN PIK3CA
TN_PIK3CA = ggplot(TN_PIK3CA, aes(x = Position, y = ratio)) +
  geom_point() +
  scale_y_continuous(expand = c(0, 0), 
                     limits = c(-1.5, 1.5)) +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  geom_hline(yintercept = seq(-1.5, 1.5, by = 0.5),
             linetype = 'dashed', size = 0.2) +
  theme(axis.line.y = element_line(size = 0.7, color = 'black'),
        axis.ticks.y = element_line(size = 0.7, lineend = 'round', color = 'black'), 
        axis.ticks.length.y = unit(0.2, 'cm'),
        panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 14, color = 'black'),
        axis.title.x = element_text(size = 14, color = 'black'),
        axis.text.x = element_blank(),
        plot.title = element_text(size = 16, color = 'black')) +
  labs(x = '', y = 'Log T/N ratio', title = 'PIK3CA')
  
TN_PIK3CA

TN_PIK3CA / DryClean_PIK3CA / Facets_PIK3CA 
  
ggsave_golden(filename = 'Figures/PIK3CA_comparision.pdf', plot = (TN_PIK3CA / DryClean_PIK3CA / Facets_PIK3CA), width = 12)

