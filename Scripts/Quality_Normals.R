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




##' replacement pattern of MODE1 DryClean
replacement = DryClean_snps[which(DryClean_snps$name == 'P-0000584-T03-IM6_P-0000584-N01'), ]
replacement$indx = seq(1, nrow(replacement), 1)
replacement$replace = as.character(as.numeric(replacement$replace))

#' chromosome borders:
co_out = c()
freq = data.frame()

for(i in unique(replacement$chrom)){
  da = replacement[which(replacement$chrom == i), ]
  co = da$indx[which.max(da$maploc)]
  co_out = c(co_out, co)
  fr = data.frame(chr = i,
                  freq = table(da$replace)[[2]] / nrow(da))
  freq = rbind(freq, fr)
  
}

genome.plot = ggplot(replacement, aes(x = indx, y = cnlr, color = replace)) +
  geom_point(size = 0.3) +
  scale_color_manual(values = c('0' = 'blue',
                                '1' = 'red')) +
  geom_vline(xintercept = co_out, color = 'grey35', size = 0.3, linetype = 'dashed') +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(-5, 5)) +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "top") +
  labs(x = 'genome', y = 'Log estimates', color = 'Replacement')

freq$chr = factor(freq$chr, levels = seq(1, 23, 1))

freq.plot = ggplot(freq, aes(x = chr, y = freq)) +
  geom_bar(stat = 'identity', color = 'black', fill = 'black') +
  geom_hline(yintercept = seq(0, 1, 0.25), linetype = 'dashed', size = 0.3, color = 'grey35') +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 1)) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank()) +
  
  labs(x = 'chromosome', y = 'substitution rate')

freq.plot

ggsave_golden(filename = 'Figures/replace_distribution.pdf', plot = (genome.plot / freq.plot), width = 12)









a = Facets_original_segments[which(Facets_original_segments$name == 'P-0016732-T02-IM6_P-0016732-N01'), ]
a = a[which(a$chrom == 17), ]
a$end - a$start








