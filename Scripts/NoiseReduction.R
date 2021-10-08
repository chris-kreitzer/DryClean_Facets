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
DataIn = readRDS('Data_out/BRCA_workingCohort_MSK.rds')
FacetsCalls = read.csv('~/Documents/MSKCC/05_IMPACT40K/Data/Signed_out/Facets_annotated.cohort.txt', sep = '\t')
sample1 = facets::readSnpMatrix('Tumor_countsFile/countsMerged____P-0000584-T03-IM6_P-0000584-N01-IM6.dat.gz')
sample1$bins = paste(sample1$Chromosome, sample1$Position, sep = ';')
sample1_cleaned = readRDS('Tumor_cleaned/countsMerged____P-0000584-T03-IM6_P-0000584-N01-IM6.dat.gz.rds_cleaned.rds')
sample1_cleaned = as.data.frame(sample1_cleaned)

RecognizedPositions = readRDS('Normal_samples/countsMerged____P-0000584-T03-IM6_P-0000584-N01-IM6.dat.gz.rds')
RecognizedPositions = as.data.frame(RecognizedPositions)
colnames(RecognizedPositions)[1] = 'chromosome'
colnames(RecognizedPositions)[2] = 'start'
RecognizedPositions$bins = paste(RecognizedPositions$chromosome, RecognizedPositions$start, sep = ';')


## Processing
TN_sample1 = sample1[which(sample1$bins %in% RecognizedPositions$bins), ]
TN_sample1$TN_ratio = log(TN_sample1$TUM.DP / TN_sample1$NOR.DP)


## Display TN ratios across chr17q (chr17: 27229517 - 48441764); including ERBB2
TN_sample1_chr17q = TN_sample1[which(TN_sample1$Chromosome == 17 & 
                                       TN_sample1$Position > 27229517 &
                                       TN_sample1$Position < 48441764), ]
TN_sample1_chr17q$indx = seq(1, nrow(TN_sample1_chr17q), 1)

ggplot(TN_sample1_chr17q, aes(x = indx, y = TN_ratio)) +
  geom_point() +
  theme(aspect.ratio = 0.5)





















## Sample2
sample2 = facets::readSnpMatrix('Tumor_countsFile/countsMerged____P-0003195-T02-IM6_P-0003195-N01-IM6.dat.gz')
sample2$bins = paste(sample2$Chromosome, sample2$Position, sep = ';')
TN_sample2 = sample2[which(sample2$bins %in% RecognizedPositions$bins), ]
TN_sample2$TN_ratio = log(TN_sample2$TUM.DP / TN_sample2$NOR.DP)
