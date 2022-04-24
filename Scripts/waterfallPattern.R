##-----------------------------------------------------------------------------
## Check WATERFALL pattern:
##-----------------------------------------------------------------------------
##
## start: 04/22/2022
## chris-kreitzer


clean()
gc()
source('~/Documents/GitHub/MSKCC/Scripts/hg19.R')
library(facets)
library(pctGCdata)
library(facetsSuite)
library(patchwork)

## P-0018672-T01-IM6
#' original
countmatrix = facets::readSnpMatrix('~/Desktop/tmp/countsMerged____P-0018672-T01-IM6_P-0018672-N01-IM6.dat.gz')
facets_original = facetsSuite::run_facets(read_counts = countmatrix,
                                          cval = 100,
                                          snp_nbhd = 250,
                                          genome = 'hg19', seed = 100)

original = cnlr_plot(facets_data = facets_original, genome = hg19)


#' Cleaned
cleaned = FacetsDC::run_facets_cleaned(read_counts = '~/Desktop/tmp/countsMerged____P-0018672-T01-IM6_P-0018672-N01-IM6.dat.gz', 
                                       read_cleaned = '~/Desktop/tmp/sample249.rds.rds', 
                                       MODE = 'union', 
                                       cval = 100, 
                                       seed = 100, snp_nbhd = 250)

cleaned_plot = cnlr_plot(facets_data = cleaned, genome = hg19)

original / cleaned_plot

## other way for visualization: 
plot(density(facets_original$snps$cnlr),
     xaxt = 'n', las = 2, 
     xlab = '',
     main = '',
     ylim = c(0, 0.8),
     lwd = 2,
     col = 'black')
axis(1, at = c(-4, -2, 0, 2, 4))
lines(density(cleaned$snps$cnlr),
      lwd = 2, col = 'red')
box(lwd = 2)
mtext(side = 1, text = 'CnLR', line = 1.8, cex = 1.3)
text(x = 1.9, y = 0.4, label = '(drycleaned-) CnLR\ndlrs=0.75', col = 'red', font = 2)
text(x = -1.8, y = 0.3, label = '(raw) CnLR\ndlrs=1.44', col = 'black', font = 2)

dlrs(facets_original$snps$cnlr)
dlrs(cleaned$snps$cnlr)


##-----------------------------------------------------------------------------
## Systematic evaluation of the waterfall pattern: 
## ----------------------------------------------------------------------------
MasterFile = readRDS('~/Documents/MSKCC/07_FacetsReview/DryClean/TUMOR_BRCA/tumor_table.rds')


countfiles = list.files('~/Documents/MSKCC/07_FacetsReview/DryClean/Tumor_Countfiles/', full.names = T)
cleanedfiles = list.files('~/Documents/MSKCC/07_FacetsReview/DryClean/TUMOR_DECOMPOSED/', full.names = T)

waterfall_out = data.frame()
for(i in 1:nrow(MasterFile)){
  try({
    print(MasterFile$original[i])
    count = facetsSuite::read_snp_matrix(input_file = grep(pattern = MasterFile$original[i], x = countfiles, value = T))
    facets_out = facetsSuite::run_facets(read_counts = count, cval = 100, genome = 'hg19', seed = 100)
    facets_dlrs = dlrs(x = facets_out$snps$cnlr)
    
    #' cleaned
    cleaned = FacetsDC::run_facets_cleaned(read_counts = grep(pattern = MasterFile$original[i], x = countfiles, value = T),
                                           read_cleaned = grep(pattern = paste0(MasterFile$sample[i], '.rds'), x = cleanedfiles, value = T),
                                           MODE = 'union', cval = 100, seed = 100)
    cleaned_dlrs = dlrs(x = cleaned$snps$cnlr)
    
    waterfall = data.frame(sample = MasterFile$original[i],
                          facets_dlrs = facets_dlrs,
                          dryclean_dlrs = cleaned_dlrs)
    
    waterfall_out = rbind(waterfall_out, waterfall)
  })
  
}

write.table(x = waterfall_out, file = 'Data_out/waterfall_dlrs.txt', sep = '\t', quote = F, row.names = F)


## Visualization:
waterfall = read.csv('Data_out/Hypersegmentation.txt', sep = '\t')

ggplot() +
  geom_jitter(aes(y = waterfall$waterfall[which(waterfall$algorithm == 'facets')], x = 1),
              position = position_jitter(width = 0.2, height = 0.1), size = 0.65, color = 'grey15', alpha = 0.35) +
  geom_boxplot(aes(y = waterfall$waterfall[which(waterfall$algorithm == 'facets')], x = 1), width = 0.2, outlier.shape = NA, color = '#850404') +
  geom_jitter(aes(y = waterfall$waterfall[which(waterfall$algorithm == 'dryclean')], x = 2), 
              position = position_jitter(width = 0.2, height = 0.1), size = 0.65, color = 'grey15', alpha = 0.35) +
  geom_boxplot(aes(y = waterfall$waterfall[which(waterfall$algorithm == 'dryclean')], x = 2), width = 0.2, outlier.shape = NA, color = '#850404') +
  scale_y_continuous(expand = c(0.05, 0.01), limits = c(0, 2.5)) +
  scale_x_continuous(expand = c(0.05, 0.01), labels = c('Facets\n(raw)', 'Facets\n(cleaned)'), breaks = c(1,2)) +
  theme_bw() +
  theme(aspect.ratio = 1.5,
        axis.text = element_text(size = 12, colour = 'black'),
        axis.title = element_text(size = 12, colour = 'black'),
        panel.border = element_rect(fill = NA, size = 1.3, color = 'black'),) +
  annotate('text', x = 1.5, y = 2.3, label = '(paired t-test) p<2.2e-16') +
  labs(x = '', y = 'dlrs [CnLR]')

t.test(waterfall$waterfall ~ waterfall$algorithm, paired = T)
