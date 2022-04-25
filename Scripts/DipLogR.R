##-----------------------------------------------------------------------------
## DipLogR and Ploidy
##-----------------------------------------------------------------------------
##
## start: 04/25/2022
## chris-kreitzer


clean()
gc()

source('~/Documents/GitHub/DryClean_Facets/Scripts/UtilityFunctions.R')
source('~/Documents/GitHub/MSKCC/Scripts/hg19.R')

sample_brca = readRDS('~/Documents/MSKCC/07_FacetsReview/DryClean/PON_BRCA/normal_table.rds')
sample_brca = substr(x = sample_brca$original, start = 1, stop = 17)
sample_brca = sample_brca[!duplicated(sample_brca)]

Annotation = read.csv('~/Documents/MSKCC/dmp-2021/dmp_022022/2022_02_09/msk_impact_facets_annotated.cohort.txt.gz', sep = '\t')
Annotation = Annotation[which(Annotation$tumor_sample %in% sample_brca), ]
ploidy_over_6 = Annotation[which(Annotation$ploidy >= 6), ]


## Fetch the right files
MasterFile = readRDS('~/Documents/MSKCC/07_FacetsReview/DryClean/TUMOR_BRCA/tumor_table.rds')
MasterFile = MasterFile[which(MasterFile$original %in% ploidy_over_6$tumor_sample), ]

countfiles = list.files('~/Documents/MSKCC/07_FacetsReview/DryClean/Tumor_Countfiles/', full.names = T)
cleanedfiles = list.files('~/Documents/MSKCC/07_FacetsReview/DryClean/TUMOR_DECOMPOSED/', full.names = T)


## TestRun: one example: P-0022981-T01-IM6
countmatrix = facets::readSnpMatrix(filename = grep(pattern = MasterFile$original[1], x = countfiles, value = T))
facets_original = facetsSuite::run_facets(read_counts = countmatrix,
                                          cval = 100,
                                          snp_nbhd = 250,
                                          genome = 'hg19', seed = 100)


facets_original$dipLogR
facets_original$ploidy
A = cnlr_plot(facets_original, genome = hg19) +
  theme(axis.text.y = element_text(size = 10, colour = 'black')) +
  annotate('text', x = 300000000, y = 4, label = 'dipLogR = -1.74\nploidy = 9.26', size = 6)



#' Cleaned
cleaned = FacetsDC::run_facets_cleaned(read_counts = grep(pattern = MasterFile$original[1], x = countfiles, value = T), 
                                       read_cleaned = grep(pattern = paste0(MasterFile$sample[1], '.rds.*'), x = cleanedfiles, value = T), 
                                       MODE = 'union', 
                                       cval = 100, 
                                       seed = 100, snp_nbhd = 250)


cleaned$dipLogR
cleaned$ploidy
B = cnlr_plot(facets_data = cleaned, genome = hg19) +
  theme(axis.text.y = element_text(size = 10, color = 'black')) +
  annotate('text', x = 300000000, y = 4, label = 'dipLogR = -0.29\nploidy = 2.81', size = 6)

A / B



ggplot(Annotation, aes(x = ploidy, y = dipLogR)) +
  geom_point(size = 0.85) +
  geom_smooth(method = 'lm', show.legend = T, inherit.aes = T) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(-1.5, 0.5)) +
  theme_bw() +
  theme(axis.text = element_text(size = 14, color = 'black'),
        panel.border = element_rect(fill = NA, color = 'black', size = 1.45, linetype = 7),
        axis.title = element_text(size = 14, color = 'black'),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        aspect.ratio = 1) +
  annotate('text', x = 6.2, y = 0.3, label = 'Pearson cor = -0.9', size = 6) +
  labs(x = 'Ploidy', y = 'dipLogR')




plot(Annotation$dipLogR ~ Annotation$ploidy)
abline(lm(Annotation$dipLogR ~ Annotation$ploidy))
cor.test(Annotation$dipLogR, Annotation$ploidy)





