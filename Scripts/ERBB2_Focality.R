## Focality of Events; look at ERBB2
##
## compare Normal Sample (quality); raw SNPs; Facets Output and summary stats 
## 
## 10/27/2021
## chris-kreitzer


set.seed(100)
source('~/Documents/GitHub/DryClean_Facets/Scripts/UtilityFunctions.R')
library(patchwork)
samples = unique(Facets_original_QC$name)


focality = function(sample, snps, icn, info, chrom, start, end, gene.start, gene.end){
  dat.info = info[which(info$sample == sample), ]
  dat.normal = snps[which(snps$name == sample), ]
  dat.focal = dat.normal[which(dat.normal$chrom == chrom & dat.normal$maploc >= start & dat.normal$maploc <= end), ]
  dat.focal$uniformity = log(dat.focal$rCountN / mean(dat.focal$rCountN))
  dat.focal$indx = seq(1, nrow(dat.focal), 1)
  
  #' general coverage throughout the whole genome
  normal.coverage = ggplot(dat.normal, aes(x = rCountN)) +
    geom_histogram(bins = 30, color = 'black', fill = 'white') +
    scale_y_continuous(expand = c(0, 0)) +
    theme_bw(base_size = 14) +
    theme(panel.grid = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(x = '', title = sample)
  
  
  #' plot normal coverage (uniformity)
  normal.focal = ggplot(dat.focal, aes(x = indx, y = uniformity)) +
    geom_jitter() +
    scale_y_continuous(expand = c(0, 0),
                       limits = c(-2.5, 2.5)) +
    geom_hline(yintercept = seq(-2.5, 2.5, 0.5), 
               col = 'grey35', 
               size = 0.2, linetype = 'dashed') +
    theme_bw(base_size = 14) +
    theme(axis.text.x = element_blank(),
          panel.background = element_rect(fill = 'grey95', colour = 'black', size = 0.8),
          axis.text.y = element_text(size = 12),
          panel.grid = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(x = '', y = 'log(N/mean(N))', title = paste0('dlrs: ', round(dlrs(dat.focal$uniformity), 3)))
    
  
  #' focality at ERBB2
  ERBB2_cnlr = ggplot(dat.focal, aes(x = indx, y = cnlr)) +
    geom_jitter() +
    scale_y_continuous(expand = c(0, 0),
                       limits = c(-2.5, 2.5)) +
    geom_hline(yintercept = seq(-2.5, 2.5, 0.5), 
               col = 'grey35', 
               size = 0.2, linetype = 'dashed') +
    geom_hline(yintercept = 0, size = 1, color = 'black') +
    theme_bw(base_size = 14) +
    theme(axis.text.x = element_blank(),
          panel.background = element_rect(fill = 'grey95', colour = 'black', size = 0.8),
          axis.text.y = element_text(size = 12),
          panel.grid = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(x = '', y = 'CnLR', title = paste0('dlrs: ', round(dlrs(dat.focal$cnlr), 3), ' @chr17_q'))
  
  
  #' extract important measures
  dat.ERBB2 = icn[which(icn$name == sample), ]
  seg = dat.ERBB2[which(dat.ERBB2$chrom == chrom),, drop = F]
  ERBB2 = seg[which(seg$start <= start & seg$end <= end & seg$end >= start | 
              seg$start <= end & seg$start >= start & seg$end >= end | 
              start <= seg$start & end >= seg$end),, drop = F]
  
  #' plot segmentation at ERBB2
  ii = cumsum(c(0, ERBB2$start))
  ii = ii / 1000000
  kk = length(ii)
  
  ERBB2.plot = ggplot(ERBB2) +
    geom_segment(aes(x = ii[-kk] + 1, xend = ii[-1], 
                     y = tcn.em, yend = tcn.em), size = 2, 
                 color = '#589E77') +

    geom_segment(aes(x = ii[-kk] + 1, xend = ii[-1], 
                     y = lcn.em, yend = lcn.em), size = 2, 
                 color = 'black') +
    
    scale_y_continuous(expand = c(0, 0),
                       limits = c(0, 7),
                       breaks = c(0, 1, 2, 3, 4, 5, 6,7)) +
    
    geom_hline(yintercept = seq(0, 7, 1), 
               col = 'grey35', 
               size = 0.2, linetype = 'dashed') +
    
    theme_bw(base_size = 14) +
    theme(axis.text.x = element_blank(),
          panel.grid = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(x = '', y = 'Integer Copy Number', title = paste0('ICN @ ', chrom, '_q'))
  
  #' Info on sample:
  sample = substr(sample, start = 1, stop = 17)
  
  ERBB2.point = seg[which(seg$start <= gene.start & seg$end <= gene.end & seg$end >= gene.start | 
                      seg$start <= gene.end & seg$end >= gene.end | 
                      gene.start <= seg$start & gene.end >= seg$end),, drop = F]
  
  tcn = ERBB2.point$tcn.em
  lcn = ERBB2.point$lcn.em
  length = round((ERBB2.point$end - ERBB2.point$start) / 1000000, 3)
  num.mark = ERBB2.point$num.mark
  nhet = ERBB2.point$nhet
  purity = round(dat.info$purity, 3)
  ploidy = round(dat.info$ploidy, 3)
  
  text_plot = ggplot() +
    annotate("text",
             x = 1,
             y = 1,
             size = 6,
             label = paste0(sample, '\n', 'ERBB2\n', 'CN = ', tcn, ':', lcn, '\n length: ', 
                            length, 'Mb\n', 'markers: ', num.mark, '\n', 'het SNPs: ', nhet, '\n',
                            'Purity: ', purity, '\n',
                            'Ploidy: ', ploidy)) + 
    theme_void()
    
    
  facets_all = list(normal.coverage, normal.focal, ERBB2_cnlr, ERBB2.plot, text_plot)
  return(facets_all)

}

gene.start = 37844167
gene.end = 37886679
start = 15908522
end = 81195210

#' applying a for loop to every sample; 
for(i in unique(samples)){
  facets = focality(sample = i, 
               snps = Facets_original_snps,
               icn = Facets_original_segments, 
               info = Facets_original_summary,
               chrom = 17, 
               start = start, 
               end = end,
               gene.start = gene.start,
               gene.end = gene.end)
  
  dryclean = focality(sample = i, 
                      snps = DryClean_snps,
                      icn = DryClean_segments, 
                      info = DryClean_summary,
                      chrom = 17, 
                      start = start, 
                      end = end,
                      gene.start = gene.start,
                      gene.end = gene.end)
  
  dryclean_full = focality(sample = i, 
                           snps = DryClean_snps_full,
                           icn = DryClean_segments_full, 
                           info = DryClean_summary_full,
                           chrom = 17, 
                           start = start, 
                           end = end,
                           gene.start = gene.start,
                           gene.end = gene.end)
  
  ggsave_golden(filename = paste0('Figures/Focal_ERBB2/', i, '.pdf'), 
                plot = plot_grid(plotlist = facets, nrow = 1) / 
                  plot_grid(plotlist = dryclean, nrow = 1) / 
                  plot_grid(plotlist = dryclean_full, nrow = 1), width = 14)
  
}


