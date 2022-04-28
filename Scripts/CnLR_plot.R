##-----------------------------------------------------------------------------
## Customized CnLR plot: based on FacetsSuite
##
## start: 04/22/2022
## 
## chris-kreitzer


source('~/Documents/GitHub/MSKCC/Scripts/hg19.R')


# cnlr_plot(facets_data = facets_data, genome = hg19)

#' function
cnlr_plot = function(facets_data,
                     colors = c('#5c9e75', '#272425'),
                     genome = c('hg19', 'hg18', 'hg38')) {
  
  genome = genome
  
  snps = facets_data$snps
  segs = facets_data$segs
  dipLogR = facets_data$dipLogR
  
  snps = get_cum_chr_maploc(snps, genome)
  mid = snps$mid[names(snps$mid) %in% snps$snps$chrom]
  centromeres = snps$centromeres
  snps = snps$snps
  
  snps$cnlr_median = rep(segs$cnlr.median, segs$num.mark)
  
  starts = cumsum(c(1, segs$num.mark))[seq_along(segs$num.mark)]
  ends = cumsum(c(segs$num.mark))
  my_starts = snps[starts, c('chr_maploc', 'cnlr_median')]
  my_ends = snps[ends, c('chr_maploc', 'cnlr_median')]
  
  ymin = floor(min(segs$cnlr.median, na.rm = T))
  ymax = ceiling(max(segs$cnlr.median, na.rm = T))
  if (ymin > -5) ymin = -5
  if (ymax < 5) ymax = 5
  
  
  pt_cols = colors[c(snps$chrom %% 2) + 1]
  
  # plot
  cnlr = ggplot(snps) +
    geom_point(aes(y = cnlr, x = chr_maploc), pch = 19, col = pt_cols, size = .4) +
    scale_x_continuous(breaks = mid, labels = names(mid), expand = c(0, 0)) +
    scale_y_continuous(breaks = scales::pretty_breaks(), limits = c(ymin, ymax), expand = c(0,0)) +
    geom_hline(yintercept = dipLogR, color = 'sandybrown', size = .8) +
    geom_segment(data = segs, aes(x = my_starts$chr_maploc, xend = my_ends$chr_maploc,
                                  y = my_starts$cnlr_median, yend = my_ends$cnlr_median),
                 col = 'red3', size = 1, lineend = 'butt') +
    labs(x = NULL, y = 'Copy number log ratio') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 0, size = 10, color = 'black'),
          axis.text.y = element_blank(),
          text = element_text(size = 12),
          axis.ticks.y = element_blank(),
          panel.grid.minor = element_line(colour = 'grey', size = 0.2),
          panel.grid.major = element_blank(),
          panel.border = element_rect(fill = NA, size = 1.5))
  return(cnlr)
  
}



get_cum_chr_maploc = function(snps, genome = hg19) {
  
  genome = genome
  
  cum_chrom_lengths = cumsum(as.numeric(genome$size))
  mid = cum_chrom_lengths - (genome$size / 2)
  names(mid) = seq_len(nrow(genome))
  centromeres = genome$centromere + c(0, cum_chrom_lengths[-length(cum_chrom_lengths)])
  
  snps$chr_maploc = snps$maploc + c(0, cum_chrom_lengths)[snps$chrom]
  
  list(snps = snps, mid = mid, centromeres = centromeres)
}
