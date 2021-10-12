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
ROI = data.frame(chromosome = c(17, 8, 1),
                 start = c(27229517, 1300569, 131354907),
                 end = c(48441764, 44119298, 239198818),
                 gene = c('ERBB2', '8p_arm', '1q_arm'))


#' TN metrics;
TN_metrics = function(data){
  file_name = basename(data)
  countdata = facets::readSnpMatrix(data, err.thresh = 10, del.thresh = 10)
  countdata$probes = paste(countdata$Chromosome, countdata$Position, sep = ';')
  countdata = countdata[which(countdata$Chromosome == ROI[1, 'chromosome']), ]
  countdata = countdata[which(countdata$probes %in% RecognizedProbes$probes[which(RecognizedProbes$chromosome == ROI[1, 'chromosome'])]), ]
  countdata = countdata[which(countdata$Position > ROI[1, 'start'] &
                                countdata$Position < ROI[1, 'end']), ]
  countdata$TN_ratio = log(countdata$TUM.DP / countdata$NOR.DP)
  countdata$indx = seq(from = 1, to = nrow(countdata), by = 1)
  countdata$sample = file_name
  countdata$dispersion = dlrs(x = countdata$TN_ratio)
  countdata
}


files = list.files('Tumor_countsFile/', full.names = T, all.files = F, include.dirs = T)
files = paste0('/Users/chriskreitzer/Documents/GitHub/DryClean_Facets/', files)


x = lapply(files, function(x) TN_metrics(x))
y = data.table::rbindlist(x)

#' 










