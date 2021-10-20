## Case study, comparing P-0000584-T03-IM6
## Comparing Facets and DryClean
##

rm(list = ls())
.rs.restartR()
set.seed(99)
library(pctGCdata)
library(patchwork)
library(DNAcopy)

#' Facets:
Sample1_Facets = facets::readSnpMatrix('Tumor_countsFile/countsMerged____P-0000584-T03-IM6_P-0000584-N01-IM6.dat.gz')
preSample1_Facets = facets::preProcSample(rcmat = Sample1_Facets)
postSample1_Facets = facets::procSample(x = preSample1_Facets)
outSample1_Facets = facets::emcncf(postSample1_Facets)

Facets_chr17 = outSample1_Facets$cncf[which(outSample1_Facets$cncf$chrom == 17), ]
Facets_chr17q = Facets_chr17[which(Facets_chr17$start > 27229517 & Facets_chr17$end < 48441764), ]

#' Visualization:
data_plot_facets1 = postSample1_Facets$jointseg[which(postSample1_Facets$jointseg$chrom == 17), ]
data_plot_facets1$indx = seq(from = 1, to = nrow(data_plot_facets1), by = 1)

plot_raw_facets = ggplot(data_plot_facets1, aes(x = indx, y = cnlr)) + 
  geom_point() +
  theme(aspect.ratio = 0.5,
        axis.line.y = element_line(colour = 'black', size = 0.2),
        panel.background = element_blank()) +
  scale_y_continuous(limits = c(-2, 2)) +
  geom_hline(yintercept = seq(-2, 2, 1), size = 0.1, linetype = 'dashed') +
  labs(title = paste0('chromosome 17'), x = 'genomic coordiantes', y = 'Facets GC-corrected CnLR') +
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
  

plot_distribution_facets = ggplot(data_plot_facets1, aes(x = cnlr)) + 
  geom_histogram(aes(y = ..density..), bins = 50, colour = "black", fill = "white") +
  labs(title = paste0('median = ', round(median(data_plot_facets1$cnlr), 3), '; sd = ', 
                      round(sd(data_plot_facets1$cnlr), 3)),
       x = 'Facets GC-corrected CnLR') +
  theme(aspect.ratio = 0.5)


#' Segmentation pattern
ii = cumsum(c(0, Facets_chr17$num.mark))
kk = length(ii)

plot_segmentation_facets = ggplot(Facets_chr17) +
  geom_segment(aes(x = ii[-kk] + 1, xend = ii[-1], y = tcn.em, yend = tcn.em),
               size = 2, color = '#589E77') +
  geom_segment(aes(x = ii[-kk] + 1, xend = ii[-1], 
                   y = lcn.em,
                   yend = lcn.em),
               size = 1.2, color = 'black') +
  annotate(
    "text", 
    label = "italic(ERBB2)", parse = T,
    x = ii[6], 
    y = 3.5, 
    size = 5, 
    colour = "red"
  ) +
  scale_y_continuous(limits = c(0, 6), breaks = c(0,1,2,3,4,5,6)) +
  theme(aspect.ratio = 0.5,
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_rect(fill = NA)) +
  labs(title = 'Segmentation pattern; chromosome 17', x = '', y = 'Integer Copy Number')

plot_segmentation_facets  

#' comprehensive output
patch1 = plot_raw_facets + plot_distribution_facets + plot_segmentation_facets
patch1 + plot_annotation(
  title = 'Facets-(default) run on P-0000584-T03-IM6',
  subtitle = 'concentrating exclusively on Chromosome 17 (ERBB2)',
  theme = theme(plot.title = element_text(size = 18)) & 
  theme(text = element_text('mono'))
)




#' concentrate on DryCleans output
sample1_DryClean = as.data.frame(readRDS('Tumor_cleaned/P-0000584-T03-IM6_P-0000584-N01-IM6_drycleaned.rds'))
sample1_DryClean_chr17 = sample1_DryClean[which(sample1_DryClean$seqnames == 17), ]

sample1_DryClean_chr17$indx = seq(from = 1, to = nrow(sample1_DryClean_chr17), by = 1)

plot_raw_dryclean = ggplot(sample1_DryClean_chr17, aes(x = indx, y = foreground.log)) + 
  geom_point() +
  theme(aspect.ratio = 0.5,
        axis.line.y = element_line(colour = 'black', size = 0.2),
        panel.background = element_blank()) +
  scale_y_continuous(limits = c(-2, 2)) +
  geom_hline(yintercept = seq(-2, 2, 1), size = 0.1, linetype = 'dashed') +
  labs(title = paste0('chromosome 17'), x = 'genomic coordiantes', y = 'DryClean foreground.log') +
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())


plot_distribution_dryclean = ggplot(sample1_DryClean_chr17, aes(x = foreground.log)) + 
  geom_histogram(aes(y = ..density..), bins = 50, colour = "black", fill = "white") +
  labs(title = paste0('median = ', round(median(sample1_DryClean_chr17$foreground.log), 3), '; sd = ', 
                      round(sd(sample1_DryClean_chr17$foreground.log), 3)),
       x = 'DryClean foreground.log') +
  theme(aspect.ratio = 0.5)


## Segmentation with DNACopy:
CBS_DryClean = sample1_DryClean[, c('seqnames', 'start', 'foreground.log')]
CNA.object = CNA(genomdat = cbind(CBS_DryClean$foreground.log),
                 chrom = CBS_DryClean$seqnames,
                 maploc = CBS_DryClean$start,
                 data.type = 'logratio',
                 sampleid = 'Sample1')
CNA.object_smoothed = smooth.CNA(x = CNA.object)
DryClean_segmented = DNAcopy::segment(x = CNA.object_smoothed,
                                      min.width = 5, undo.splits = 'sdundo', undo.SD = 3)


DryClean_segmented_chr17 = DryClean_segmented$output[which(DryClean_segmented$output$chrom == 17), ]
oo = cumsum(c(0, DryClean_segmented_chr17$num.mark))
jj = length(oo)

plot_segmentation_dryclean = ggplot(DryClean_segmented_chr17) +
  geom_segment(aes(x = oo[-jj] + 1,
                   xend = oo[-1],
                   y = seg.mean,
                   yend = seg.mean),
               size = 2, 
               color = '#589E77') +
  scale_y_continuous(limits = c(-1.5, 1.5), breaks = c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5)) +
  theme(panel.border = element_rect(fill = NA),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(title = 'CBS of DryClean foreground-log',
       x = '', y = 'segmentation mean') +
  annotate(
    "text", 
    label = "italic(ERBB2)", parse = T,
    x = oo[25], 
    y = 0.9, 
    size = 4, 
    colour = "red"
  )
  

#' patch2:
patch2 = plot_raw_dryclean + plot_distribution_dryclean + plot_segmentation_dryclean
patch2 = patch2 + plot_annotation(
  title = 'DryClean run on P-0000584-T03-IM6 WITH CBS',
  theme = theme(plot.title = element_text(size = 18)) & 
    theme(text = element_text('mono'))
)



patch1 / patch2

#' out





