## DryClean; Creating Panel of Normal:
## 
## Here, we will start focusing on Breast adenocarcinoma samples
## We will be using IMPACT gene panel 468 for this appraoch;
## Our PON will consist of 100 normal samples
## 
## Full detail: Data queried at 07/07/2021 from cBIO
## Breast Cancer; Breast Invasive Ductal Carcinoma; Gene Panel 468 (n = 2,705 samples)
## 
## 07/07/2021
## chris kreitzer


rm(list = ls())
.rs.restartR()
setwd('~/Documents/MSKCC/07_FacetsReview/')
set.seed(111)


## Libraries and Dependencies
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = TRUE)
# devtools::install_github('mskilab/dryclean')
# devtools::install_github('mskilab/gUtils')
# devtools::install_github('mskilab/fragCounter')
# devtools::install_github("mskilab/bamUtils")
# devtools::install_github('mskilab/fragCounter')   
# BiocManager::install('S4Vectors')
# BiocManager::install('GenomicAlignments')
# library(S4Vectors)
library(gUtils)
library(dryclean)


## Input data:
BRCA = read.csv('Data4Analysis/Breast_clinicalData_07.07.21.tsv', sep = '\t')
Paths = read.csv('~/Documents/MSKCC/05_IMPACT40K/Data/Signed_out/Facets_annotated.cohort.txt', sep = '\t')
BRCA_PON_list = readRDS('DataProcessed/BRCA_PON_list.rds')


## Data wrangling and processing
#' create a Panel of Normal
random.Normals = BRCA[sample(nrow(BRCA), size = 100, replace = F), 'Sample.ID']
paths.Normals = Paths[which(Paths$tumor_sample %in% random.Normals), 'counts_file']
write.table(paths.Normals, file = 'DataProcessed/PON_BRCA_Paths.txt', col.names = F, row.names = F, quote = F)

#' fetch coordinates from Normal samples; 
#' this script will be submitted to LFS on the juno-cluster
# library(facets)
# PON = read.csv('/juno/home/kreitzec/DryClean/PON_BRCA_Paths.txt', sep = '\t', header = F, row.names = F)
# 
# BRCA_PON_list = list()
# for(i in unique(PON)){
#   data.in = facets::readSnpMatrix(i)
#   data.processed = data.in[which(data.in$NOR.DP >= 35), c('Chromosome', 'Position', 'NOR.DP', 'NOR.RD')]
#   data.processed$sample = substr(x = basename(i), start = 17, stop = 33)
#   BRCA_PON_list[[i]] = data.processed
# }
# 
# saveRDS(object = BRCA_PON_list, file = '/juno/home/kreitze/DryClean/BRCA_PON_list.rds')


BRCA_PON_list = readRDS('DataProcessed/BRCA_PON_list.rds')
BRCA_PON_df = do.call('rbind', BRCA_PON_list)

#' automate marker selection for proper dimensions in PON
#' Note, that n (marker-bins) x m(samples) need to be equal among all normal samples

prepare_PON = function(normal_samples, input_format = NULL, sample_threshold = NULL){
  #' if the number of postions should be increased;
  #' we can introduce a sub-function here that all the postions which are present in 90% of samples
  #' and fill the remaining positions with 0;
  #' with the mean normalization, this should not influence the overall normalization
  bins_PON = data.frame()
  sample_PON = data.frame()
  message('Welcome. I am creating an n x m matrix which is equal among all normal samples')
  
  if(!is.null(input_format)){
    message('Input is converted to data.frame')
    input_list = do.call('rbind', normal_samples)
  } else{
    message('Input is data.frame')
    input_list = normal_samples
  }
  
  #' loop through list and select common positions
  for(chromosome in c(as.character(seq(1, 22, 1)), 'X')){
    chromosome_subset = input_list[which(input_list$Chromosome == chromosome), ]
    chromosome_position = as.data.frame(table(chromosome_subset$Position))
    chromosome_position$Chromosome = chromosome
    positions_keep = data.frame(loc = chromosome_position$Var1[which(chromosome_position$Freq == 99)],
                                chromosome = chromosome)
    
    #' subset normal PON
    input_selected = input_list[which(input_list$Chromosome == chromosome & input_list$Position %in% positions_keep$loc), c('Chromosome','Position', 'NOR.DP', 'sample')]
    input_selected$bin = paste(input_selected$Chromosome, input_selected$Position, sep = ';')
    sample_PON = rbind(sample_PON, input_selected)
    bins_PON = rbind(bins_PON, positions_keep)
  }
  
  return(list(bins_PON,
              sample_PON))
}

PON_out = prepare_PON(normal_samples = BRCA_PON_df)

#' downstream modification
PON_out = as.data.frame(do.call('cbind', split(y[, c('NOR.DP')], y$sample)))
row.names(PON_out) = y$bin[which(y$sample == 'P-0028201-T01-IM6')]



#' make a cross-check if matrix is okay;
#' start with mean-normalization
#' start prepare_detergent - DryClean (and documentation)












