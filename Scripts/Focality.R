## Overall Focality comparision between the three modes
## compare segment lengths from i) Facets calls ii) partial Dryclean and iii) full DryClean
## 


focal = c(0, 10)
focal_int = c(10, 50)
inter = c(50, 200)
broad = c(200, 1000)
large = c(>1)



#' units in Kb
Facets_focal = Facets_original_segments
Facets_focal$length = (Facets_focal$end - Facets_focal$start) / 1000
Facets_focal$type = NA

Facets_focal$type = ifelse(between(Facets_focal$length, focal[1], focal[2]), 'focal',
                           ifelse(between(Facets_focal$length, focal_int[1], focal_int[2]), 'focal_int',
                                  ifelse(between(Facets_focal$length, inter[1], inter[2]), 'inter',
                                                 ifelse(between(Facets_focal$length, broad[1], broad[2]), 'broad', 'large'))))

Facets_focality = data.frame(table(Facets_focal$type))
Facets_focality$rel = (Facets_focality$Freq / sum(table(Facets_focal$type))[1]) * 100
Facets_focality$type = 'Facets'


#' partial replacement
DryClean_focal = DryClean_segments
DryClean_focal$length = (DryClean_focal$end - DryClean_focal$start) / 1000
DryClean_focal$type = NA

DryClean_focal$type = ifelse(between(DryClean_focal$length, focal[1], focal[2]), 'focal',
                           ifelse(between(DryClean_focal$length, focal_int[1], focal_int[2]), 'focal_int',
                                  ifelse(between(DryClean_focal$length, inter[1], inter[2]), 'inter',
                                         ifelse(between(DryClean_focal$length, broad[1], broad[2]), 'broad', 'large'))))
sum(table(DryClean_focal$type))
Dryclean_focality = data.frame(table(DryClean_focal$type))
Dryclean_focality$rel = (Dryclean_focality$Freq / sum(table(DryClean_focal$type))[1]) * 100
Dryclean_focality$type = 'DryClean[:partial]'


#' full replacement
DryClean_focal_full = DryClean_segments_full
DryClean_focal_full$length = (DryClean_focal_full$end - DryClean_focal_full$start) / 1000
DryClean_focal_full$type = NA

DryClean_focal_full$type = ifelse(between(DryClean_focal_full$length, focal[1], focal[2]), 'focal',
                             ifelse(between(DryClean_focal_full$length, focal_int[1], focal_int[2]), 'focal_int',
                                    ifelse(between(DryClean_focal_full$length, inter[1], inter[2]), 'inter',
                                           ifelse(between(DryClean_focal_full$length, broad[1], broad[2]), 'broad', 'large'))))

Dryclean_full_focality = data.frame(table(DryClean_focal_full$type))
Dryclean_full_focality$rel = (Dryclean_full_focality$Freq / sum(table(DryClean_focal_full$type))[1]) * 100
Dryclean_full_focality$type = 'DryClean[:full]'


#' merge the data frames and make a bar_chart
Focality_Assessment = rbind(Facets_focality, Dryclean_focality, Dryclean_full_focality)
Focality_Assessment$Var1 = factor(Focality_Assessment$Var1, levels = rev(c('large', 'broad', 'inter', 'focal_int', 'focal')))
Focality_Assessment$type = factor(Focality_Assessment$type, levels = c('Facets', 'DryClean[:partial]', 'DryClean[:full]'))

focal_full = ggplot(Focality_Assessment, aes(x = type, y = rel, fill = Var1)) +
  geom_bar(stat = 'identity', position = 'stack') +
  scale_fill_manual(values = c('large' = '#0F4784',
                               'broad' = '#F04722',
                               'inter' = '#FDD31E',
                               'focal_int' = '#579E43',
                               'focal' = '#7D1123'),
                    name = '',
                    labels = c('> 1Mb', '200-1000 kb', '50-200 kb', '10-50 kb', '<10 kb')) +
  scale_y_continuous(expand = c(0,0)) +
  #scale_x_discrete(expand = c(0.1, 0.1)) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = 'black', size = 16),
        aspect.ratio = 1,
        legend.position = 'top') +
  labs(x = '', y = 'Percent')



focal_partial = ggplot(Focality_Assessment[which(Focality_Assessment$Var1 %in% c('broad', 'inter', 'focal_int', 'focal')), ], 
                       aes(x = type, y = rel, fill = Var1)) +
  geom_bar(stat = 'identity', position = 'stack') +
  scale_fill_manual(values = c('broad' = '#F04722',
                               'inter' = '#FDD31E',
                               'focal_int' = '#579E43',
                               'focal' = '#7D1123'),
                    name = '',
                    labels = c('200-1000 kb', '50-200 kb', '10-50 kb', '<10 kb')) +
  scale_y_continuous(expand = c(0,0)) +
  #scale_x_discrete(expand = c(0.1, 0.1)) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = 'black', size = 16),
        aspect.ratio = 1,
        legend.position = 'top') +
  labs(x = '', y = 'Percent')

ggsave_golden(filename = 'Figures/Focality_comparison.pdf', plot = (focal_full + focal_partial), width = 14)










