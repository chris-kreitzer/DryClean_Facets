## Extract pivotal information from Facets calls; and compare 
## pre- post DryClean outputs; find marker postions where DryClean makes an impact
## 
## 10/06/2021: original script
## chris-kreitzer
## 


rm(list = ls())
.rs.restartR()
setwd('~/Documents/GitHub/DryClean_Facets/')
set.seed(12)


## Take the T/N ratio of the same positions as used for DryClean and look into
## the distribution (deviation; noise reduction)


## Libraries and Input;
## Starting with 1 sample (P-0000584-T03-IM6)
library(patchwork)
library(pctGCdata)

DataIn = readRDS('Data_out/BRCA_workingCohort_MSK.rds')
FacetsCalls = read.csv('~/Documents/MSKCC/05_IMPACT40K/Data/Signed_out/Facets_annotated.cohort.txt', sep = '\t')
sample1 = facets::readSnpMatrix('Tumor_countsFile/countsMerged____P-0000584-T03-IM6_P-0000584-N01-IM6.dat.gz')
sample1$bins = paste(sample1$Chromosome, sample1$Position, sep = ';')
sample1_cleaned = readRDS('Tumor_cleaned/countsMerged____P-0000584-T03-IM6_P-0000584-N01-IM6.dat.gz.rds_cleaned.rds')
sample1_cleaned = as.data.frame(sample1_cleaned)
colnames(sample1_cleaned)[1] = 'chromosome'
colnames(sample1_cleaned)[2] = 'start'

#' bins to keep:
RecognizedPositions = readRDS('Normal_samples/countsMerged____P-0000584-T03-IM6_P-0000584-N01-IM6.dat.gz.rds')
RecognizedPositions = as.data.frame(RecognizedPositions)
colnames(RecognizedPositions)[1] = 'chromosome'
colnames(RecognizedPositions)[2] = 'start'
RecognizedPositions$bins = paste(RecognizedPositions$chromosome, RecognizedPositions$start, sep = ';')


## Processing
TN_sample1 = sample1[which(sample1$bins %in% RecognizedPositions$bins), ]
TN_sample1$TN_ratio = log(TN_sample1$TUM.DP / TN_sample1$NOR.DP)

#' Facets default processing of samples
Facets_sample1 = facets::preProcSample(sample1, gbuild = 'hg19', snp.nbhd = 0)
Facets_sample1_chr17q = Facets_sample1$jointseg[which(Facets_sample1$jointseg$chrom == 17), ]
Facets_sample1_chr17q = Facets_sample1_chr17q[which(Facets_sample1_chr17q$maploc %in% TN_sample1_chr17q$Position), ]
Facets_sample1_chr17q$indx = seq(1, nrow(Facets_sample1_chr17q), 1)

## Display TN ratios across chr17q (chr17: 27229517 - 48441764); including ERBB2
TN_sample1_chr17q = TN_sample1[which(TN_sample1$Chromosome == 17 & 
                                       TN_sample1$Position > 27229517 &
                                       TN_sample1$Position < 48441764), ]
TN_sample1_chr17q$indx = seq(1, nrow(TN_sample1_chr17q), 1)
sample1_cleaned_chr17q = sample1_cleaned[which(sample1_cleaned$chromosome == 17 &
                                                 sample1_cleaned$start > 27229517 &
                                                 sample1_cleaned$start < 48441764), ]
sample1_cleaned_chr17q$indx = seq(1, nrow(sample1_cleaned_chr17q), 1)

#' mark one gene
ERBB2_coordinates = sample1_cleaned_chr17q$indx[which(sample1_cleaned_chr17q$start >= 37842338 & 
                                                        sample1_cleaned_chr17q$start <= 37886915)]

#' make a visualization
TN_raw = ggplot(TN_sample1_chr17q, aes(x = indx, y = TN_ratio)) +
  geom_point() +
  theme(aspect.ratio = 0.5,
        axis.line.y = element_line(colour = 'black', size = 0.2),
        panel.background = element_blank()) +
  scale_y_continuous(limits = c(-2, 2)) +
  geom_hline(yintercept = seq(-2, 2, 1), size = 0.1, linetype = 'dashed') +
  geom_segment(aes(x = ERBB2_coordinates[1], xend = ERBB2_coordinates[length(ERBB2_coordinates)],
                   y = 1.1, yend = 1.1), color = 'red', size = 1.2) +
  geom_text(x = ERBB2_coordinates[1], y = 1.15, label = "ERBB2", vjust = 'middle')


TN_corrected = ggplot(Facets_sample1_chr17q, aes(x = indx, y = cnlr)) +
  geom_point() +
  theme(aspect.ratio = 0.5,
        axis.line.y = element_line(colour = 'black', size = 0.2),
        panel.background = element_blank()) +
  scale_y_continuous(limits = c(-2, 2)) +
  geom_hline(yintercept = seq(-2, 2, 1), size = 0.1, linetype = 'dashed') +
  geom_segment(aes(x = ERBB2_coordinates[1], xend = ERBB2_coordinates[length(ERBB2_coordinates)],
                   y = 1.1, yend = 1.1), color = 'red', size = 1.2)

TN_cleaned = ggplot(sample1_cleaned_chr17q, aes(x = indx, y = foreground.log)) +
  geom_point() +
  theme(aspect.ratio = 0.5,
        axis.line.y = element_line(colour = 'black', size = 0.2),
        panel.background = element_blank()) +
  scale_y_continuous(limits = c(-2, 2)) +
  geom_hline(yintercept = seq(-2, 2, 1), size = 0.1, linetype = 'dashed') +
  geom_segment(aes(x = ERBB2_coordinates[1], xend = ERBB2_coordinates[length(ERBB2_coordinates)],
                   y = 1.1, yend = 1.1), color = 'red', size = 1.2)



TN_raw / TN_corrected / TN_cleaned

dlrs(x = TN_sample1_chr17q$TN_ratio)
dlrs(x = Facets_sample1_chr17q$cnlr)
dlrs(x = sample1_cleaned_chr17q$foreground.log)



cor.test(TN_sample1_chr17q$TN_ratio, sample1_cleaned_chr17q$foreground.log)



## Sample2
sample2 = facets::readSnpMatrix('Tumor_countsFile/countsMerged____P-0003195-T02-IM6_P-0003195-N01-IM6.dat.gz')
sample2$bins = paste(sample2$Chromosome, sample2$Position, sep = ';')
TN_sample2 = sample2[which(sample2$bins %in% RecognizedPositions$bins), ]
TN_sample2$TN_ratio = log(TN_sample2$TUM.DP / TN_sample2$NOR.DP)
