#' Noise detection; population wise (selected cohort)
#' Discriminate between high and low purity samples
#' 
#' n = 65 selected samples
#' 
#' 10/11/2021: adapt ideas from single run and apply to all samples
#' chris-kreitzer
#' 

rm(list = ls())
.rs.restartR()
set.seed(99)
setwd('~/Documents/GitHub/DryClean_Facets/')
source('Scripts/UtilityFunctions.R')

## Libraries and Input
DataIn = readRDS('Data_out/BRCA_workingCohort_MSK.rds')

#' probes to investigate
RecognizedProbes = as.data.frame(readRDS('Normal_samples/countsMerged____P-0000584-T03-IM6_P-0000584-N01-IM6.dat.gz.rds'))
colnames(RecognizedProbes)[1] = 'chromosome'
colnames(RecognizedProbes)[2] = 'start'
RecognizedProbes$probes = paste(RecognizedProbes$chromosome, RecognizedProbes$start, sep = ';')

#' ROI
ROI = data.frame(chromosome = c(17, ))
#' TN metrics;
TN_metrics = function(data){
  countdata = facets::readSnpMatrix(data, err.thresh = 10, del.thresh = 10)
  countdata$probes = paste(countdata$Chromosome, countdata$Position)
}

sample1 = facets::readSnpMatrix('Tumor_countsFile/countsMerged____P-0000584-T03-IM6_P-0000584-N01-IM6.dat.gz')
head(sample1)
head(RecognizedProbes)
