## Normal quality: Assessing the suitability of Normals
## 


#' Data which is required:
#' Facets_original_snps - 

plot_list = list()
for(i in unique(Facets_original_snps$name)){
  dat = Facets_original_snps[which(Facets_original_snps$name == i), ]
  dat$normalized = log(dat$rCountN / mean(dat$rCountN))
  dat$indx = seq(1, nrow(dat), 1)
  hist(dat$normalized, main = i, nclass = 30, xlab = 'log(N/mean(N))')
  plot = ggplot(dat, aes(x = indx, y = normalized)) + 
    geom_point() +
    scale_y_continuous(expand = c(0, 0),
                       limits = c(-3, 3)) +
    scale_x_continuous(expand = c(0.01, 0.01)) +
    geom_hline(yintercept = seq(-2, 2, 1), linetype = 'dashed', color = 'grey35') +
    theme_base(base_size = 14)
  plot_list[[i]] = plot
}



a = Facets_original_snps[which(Facets_original_snps$name == 'P-0000584-T03-IM6_P-0000584-N01'), ]
b = DryClean_snps[which(DryClean_snps$name == 'P-0000584-T03-IM6_P-0000584-N01'), ]


