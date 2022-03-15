##-----------------------------------------------------------------------------
## investigate chromosome 17 as an example;
library(tidyverse)
normalized = vroom::vroom('~/Documents/MSKCC/07_FacetsReview/DryClean/DataProcessed/PON_normalized.txt', delim = '\t')
normalized = as.data.frame(normalized)
ii = grep(pattern = '17;.*', x = normalized$duplication)

normal_17 = normalized[ii, ]
normal_17 = separate(normal_17,
                     col = 'duplication',
                     into = c('chromosome', 'start'),
                     sep = ';', 
                     remove = F)
normal_17$start = as.numeric(as.character(normal_17$start))

#' just work with bins on chromosome 17q
normal17q = normal_17[which(normal_17$start > 24000000), ]

#' normal seq coverage on the 17q arm; over-representation
row.names(normal17q) = normal17q$duplication
normal17q[normal17q == 1] = NA
average = apply(normal17q[, 4:ncol(normal17q)], 1, FUN = function(x) mean(x, na.rm = T))
average.df = data.frame(bin = names(average),
                        value = average)

#' which bins are over-represented >> 99% quantile
over = names(average)[which(average > quantile(average, probs = c(0.99)))]
average.df$outlier = NA
for(i in 1:nrow(average.df)){
  if(average.df$bin[i] %in% over){
    average.df$outlier[i] = 'yes'
  }
}

average.df = separate(average.df,
                      'bin', 
                      into = c('chromosome', 'start'),
                      sep = ';', 
                      remove = F)

average.df$indx = seq(1, nrow(average.df), 1)
par(mfrow = c(6,1),
    mar = c(0.2,4,0.2,1))
plot(average.df$indx, 
     average.df$value, 
     pch = 20, 
     cex = 0.35, 
     adj = 0,
     ylab = '',
     xlab = '',
     xaxt = 'n', 
     las = 2)
title(ylab = "Population: NORMAL seq. depth", line = 1.8) 
box(lwd = 2)
abline(h = quantile(average, probs = c(0.90, 0.99)), lty = 'dashed', col = 'red', lwd = 1)
points(average.df$indx[!is.na(average.df$outlier)], average.df$value[!is.na(average.df$outlier)], 
       col = 'red', cex = 0.40, pch = 20)


#' extract regions where we have elevated seq. representation
average.df$annotation = NA
average.df$annotation[which(average.df$bin == '17;29095900')] = 'SUZ12P1_exon5'
average.df$annotation[which(average.df$bin == '17;30302550')] = 'SUZ12_exon6'
average.df$annotation[which(average.df$bin == '17;33817050')] = 'SLFN12L_random'
average.df$annotation[which(average.df$bin == '17;37872647')] = 'ERBB2_exon13'
average.df$annotation[which(average.df$bin == '17;40354350')] = 'STAT5B_exon18'
average.df$annotation[which(average.df$bin == '17;40354750')] = 'STAT5B_exon17'
average.df$annotation[which(average.df$bin == '17;40364000')] = 'STAT5B_exon13'
average.df$annotation[which(average.df$bin == '17;40370200')] = 'STAT5B_exon9'
average.df$annotation[which(average.df$bin == '17;40371848')] = 'STAT5B_exon6'
average.df$annotation[which(average.df$bin == '17;40379585')] = 'STAT5B_exon3'
average.df$annotation[which(average.df$bin == '17;40441938')] = 'STAT5A_exon3'
average.df$annotation[which(average.df$bin == '17;40451750')] = 'STAT5A_exon6'
average.df$annotation[which(average.df$bin == '17;40453350')] = 'STAT5A_exon8'
average.df$annotation[which(average.df$bin == '17;40458250')] = 'STAT5A_exon13'
average.df$annotation[which(average.df$bin == '17;40461100')] = 'STAT5A_exon17'
average.df$annotation[which(average.df$bin == '17;40461450')] = 'STAT5A_exon18'
average.df$annotation[which(average.df$bin == '17;41231250')] = 'BRCA1_exon13'
average.df$annotation[which(average.df$bin == '17;41245488')] = 'BRCA1_exon10'
average.df$annotation[which(average.df$bin == '17;43345050')] = 'MAP3K14_AS1'
average.df$annotation[which(average.df$bin == '17;43364100')] = 'MAP3K14_exon6'
average.df$annotation[which(average.df$bin == '17;46805550')] = 'HOXB13_exon1'
average.df$annotation[which(average.df$bin == '17;55752391')] = 'MSI2_exon12'
average.df$annotation[which(average.df$bin == '17;70120278')] = 'SOX9_exon3'
average.df$annotation[which(average.df$bin == '17;73774679')] = 'H3-3B_exon4'

for(i in 1:nrow(average.df)){
  if(!is.na(average.df$annotation[i])){
    text(x = average.df$indx[i],
         y = average.df$value[i],
         label = average.df$annotation[i],
         pos = 4,
         cex = 0.50,
         offset = 0.1)
  }
}


##-----------------------------------------------------------------------------
## The higher the seq. coverage in the normal, the more likelier we see a diploid
## or a loss in this particular region (conventional algorithms)

#' normal coverage
samp = normalized[, c(1, grep(pattern = 'P-0024637-T01-IM6_P-0024637-N01-IM6.dat.gz;.*', x = colnames(normalized)))]
samp17q = samp[which(samp$duplication %in% normal17q$duplication),, drop = F]
samp17q$indx = seq(1, nrow(samp17q), 1)
colnames(samp17q)[2] = 'value'
plot(samp17q$indx[which(samp17q$value != 1)], 
     samp17q$value[which(samp17q$value != 1)],
     pch = 20, 
     cex = 0.35, 
     adj = 0,
     ylab = '',
     xlab = '',
     xaxt = 'n', las = 2)
title(ylab = "Sample 1: Normal seq. depth", line = 1.8) 
box(lwd = 2)
abline(h = quantile(average, probs = c(0.90, 0.99)), lty = 'dashed', col = 'red', lwd = 1)

#' tumor coverage
tumor_normalized = readRDS('Data4Analysis/TumorNormalizedAll.rds')
tumor_normalized = as.data.frame(tumor_normalized)
samp_tumor = tumor_normalized[, c(1, grep(pattern = '^P-0024637-T01-IM6.*', colnames(tumor_normalized)))]
samp_tumor17q = samp_tumor[which(samp_tumor$duplication %in% normal17q$duplication),, drop = F]
samp_tumor17q$indx = seq(1, nrow(samp_tumor17q), 1)

colnames(samp_tumor17q)[2] = 'value'
plot(samp_tumor17q$indx[which(samp_tumor17q$value != 1)], 
     samp_tumor17q$value[which(samp_tumor17q$value != 1)],
     pch = 20, 
     cex = 0.35, 
     adj = 0,
     ylab = '',
     xlab = '',
     xaxt = 'n', las = 2)
title(ylab = "Sample 1: Tumor seq. depth", line = 1.8) 
box(lwd = 2)
abline(h = quantile(average, probs = c(0.90, 0.99)), lty = 'dashed', col = 'red', lwd = 1)


#' decomposed tumor foreground - true signals:
tumor_decomposed = readRDS('Tumor_Decomposed/sample464.rds.rds')
tumor_decomposed = as.data.frame(tumor_decomposed)
tumor_decomposed$bin = paste(tumor_decomposed$seqnames, tumor_decomposed$start, sep = ';')
tumor_decomposed17q = tumor_decomposed[which(tumor_decomposed$bin %in% normal17q$duplication),, drop = F]
tumor_decomposed17q$indx = seq(1, nrow(tumor_decomposed17q), 1)

plot(tumor_decomposed17q$indx[which(tumor_decomposed17q$input.read.counts != 1)], 
     tumor_decomposed17q$foreground[which(tumor_decomposed17q$input.read.counts != 1)],
     pch = 20, 
     cex = 0.35, 
     adj = 0,
     ylab = '',
     xlab = '',
     xaxt = 'n', las = 2)
title(ylab = "Sample 1: FOREGROUND", line = 1.8) 
box(lwd = 2)
abline(h = quantile(average, probs = c(0.90, 0.99)), lty = 'dashed', col = 'red', lwd = 1)


#' decomposed tumor background - extraction:
tumor_decomposed = readRDS('Tumor_Decomposed/sample464.rds.rds')
tumor_decomposed = as.data.frame(tumor_decomposed)
tumor_decomposed$bin = paste(tumor_decomposed$seqnames, tumor_decomposed$start, sep = ';')
tumor_decomposed17q = tumor_decomposed[which(tumor_decomposed$bin %in% normal17q$duplication),, drop = F]
tumor_decomposed17q$indx = seq(1, nrow(tumor_decomposed17q), 1)

plot(tumor_decomposed17q$indx[which(tumor_decomposed17q$input.read.counts != 1)], 
     tumor_decomposed17q$background[which(tumor_decomposed17q$input.read.counts != 1)],
     pch = 20, 
     cex = 0.35, 
     adj = 0,
     ylab = '',
     xlab = '',
     xaxt = 'n', las = 2)
title(ylab = "Sample 1: BACKGROUND", line = 1.8) 
box(lwd = 2)
abline(h = quantile(average, probs = c(0.90, 0.99)), lty = 'dashed', col = 'red', lwd = 1)


#' Facets original calls
facets_original = readRDS('Tumor_Facets/countsMerged____P-0024637-T01-IM6_P-0024637-N01-IM6.dat.gz.rds')
segs = facets_original$facets_output$segs
segs = segs[which(segs$chrom == 17), ]

plot(0, 
     xaxt = 'n', 
     yaxt = 'n', 
     bty = 'n', 
     pch = '', 
     ylab = '', 
     xlab = '',
     xlim = c(1, 12522),
     ylim = c(0, 6),
     tck = 0)
axis(side = 2, at = c(0, 1, 2, 3, 4, 5), labels = c(0, 1, 2, 3, 4, 5), las = 2)
box(lwd = 2)
title(ylab = 'FACETS: Integer Copy Number', line = 1.8)
segments(x0 = 1, x1 = 11271, 
         y0 = 4, y1 = 4, 
         col = '#91b1a1', lwd = 2.5)

segments(x0 = 1, x1 = 11271, 
         y0 = 2, y1 = 2, 
         col = 'grey55', lwd = 2.5)

#' second segment
segments(x0 = 11278, x1 = 12522, 
         y0 = 1, y1 = 1, 
         col = '#91b1a1', lwd = 2.5)

#' second segment
segments(x0 = 11278, x1 = 12522, 
         y0 = 0, y1 = 0, 
         col = 'grey55', lwd = 2.5)



#' Facets dryclean calls
facets_dryclean = readRDS('Facets_DryClean/sample464.rds')
segs_clean = facets_dryclean$Facets$segs
segs_clean = segs_clean[which(segs_clean$chrom == 17), ]

plot(0, 
     xaxt = 'n', 
     yaxt = 'n', 
     bty = 'n', 
     pch = '', 
     ylab = '', 
     xlab = '',
     xlim = c(1, 12522),
     ylim = c(0, 6),
     tck = 0)
axis(side = 2, at = c(0, 1, 2, 3, 4, 5), labels = c(0, 1, 2, 3, 4, 5), las = 2)
box(lwd = 2)
title(ylab = 'Dryclean: Integer Copy Number', line = 1.8)
segments(x0 = 1, x1 = which.max(tumor_decomposed17q$indx[which(tumor_decomposed17q$start <= segs_clean$end[2])]), 
         y0 = 3, y1 = 3, 
         col = '#91b1a1', lwd = 2.5)

segments(x0 = 1, x1 = which.max(tumor_decomposed17q$indx[which(tumor_decomposed17q$start <= segs_clean$end[2])]), 
         y0 = 1, y1 = 1, 
         col = 'grey55', lwd = 2.5)

#' second segment
segments(x0 = which.max(tumor_decomposed17q$indx[which(tumor_decomposed17q$start < segs_clean$end[2])]), 
         x1 = nrow(tumor_decomposed17q), 
         y0 = 2, y1 = 2, 
         col = '#91b1a1', lwd = 2.5)

#' second segment
segments(x0 = which.max(tumor_decomposed17q$indx[which(tumor_decomposed17q$start < segs_clean$end[2])]), 
         x1 = nrow(tumor_decomposed17q), 
         y0 = 0, y1 = 0, 
         col = 'grey55', lwd = 2.5)



