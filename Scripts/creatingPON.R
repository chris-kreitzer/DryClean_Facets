## DryClean; Creating Panel of Normal:
## 
## Here, we will start focusing on Breast adenocarcinoma samples
## We will be using IMPACT gene panel 468 for this appraoch;
## Our PON will consist of 1000 normal samples
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
devtools::install_github('mskilab/dryclean')
# devtools::install_github('mskilab/gUtils')
# devtools::install_github('mskilab/fragCounter')
# devtools::install_github("mskilab/bamUtils")
# devtools::install_github('mskilab/fragCounter')   
# BiocManager::install('S4Vectors')
# BiocManager::install('GenomicAlignments')
# library(S4Vectors)
library(gUtils)
library(dryclean)
library(tidyverse)
library(pbmcapply)
library(data.table)

## Input data:
BRCA = read.csv('Data4Analysis/Breast_clinicalData_07.07.21.tsv', sep = '\t')
Paths = read.csv('~/Documents/MSKCC/05_IMPACT40K/Data/Signed_out/Facets_annotated.cohort.txt', sep = '\t')
BRCA_PON_list = readRDS('DataProcessed/BRCA_PON_list.rds')


## Data wrangling and processing
#' create a Panel of Normal; n = 1,000
random.Normals = BRCA[sample(nrow(BRCA), size = 1020, replace = F), 'Sample.ID']
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

#' replace with Rbindlist
BRCA_PON_df = rbindlist(BRCA_PON_list)

#' automate marker selection for proper dimensions in PON
#' Note, that n (marker-bins) x m(samples) need to be equal among all normal samples

prepare_PON = function(normal_samples, sample_threshold = NULL){
  #' if the number of postions should be increased;
  #' we can introduce a sub-function here that all the postions which are present in 90% of samples
  #' and fill the remaining positions with 0;
  #' with the mean normalization, this should not influence the overall normalization
  
  input_list = normal_samples
  bins_PON = data.frame()
  sample_PON = data.frame()
  message('Welcome. I am creating an n x m matrix which is equal among all normal samples')
  
  #' search for samples which have duplicated entries (chromosome*postion) and discard them
  input_list$duplication = paste(input_list$Chromosome, input_list$Position, sep = ';')
  container = c()
  for(i in unique(input_list$sample)){
    if(any(duplicated(input_list$duplication[which(input_list$sample == i)]))){
      container = c(container, i)
    } 
    else next
  }
  
  message('Samples')
  #' subset input list; to remove samples with duplicated entries
  input_list = input_list[!input_list$sample %in% container, ]
  rm(container)
  message('Sample Quality Control ended')
  
  
  #' loop through list and select common positions
  n.PON = length(unique(input_list$sample))
  threshold = round(n.PON * sample_threshold)
  matrix.table = data.frame(table(input_list$duplication))
  matrix.table$Var1 = as.character(as.factor(matrix.table$Var1))
  matrix.table.keep = matrix.table[which(matrix.table$Freq >= threshold), ]
  colnames(matrix.table.keep)[1] = 'loc'
  
  #' sample-wise listing of postions
  locations.out = data.frame()
  for(patient in unique(input_list$sample)){
    print(patient)
    if(all(matrix.table.keep$loc %in% input_list$duplication[which(input_list$sample == patient)])){
      table.out = input_list[which(input_list$sample == patient & input_list$duplication %in% matrix.table.keep$loc), ]
    } else {
      table.out = input_list[which(input_list$sample == patient & input_list$duplication %in% matrix.table.keep$loc), ]
      missing = setdiff(matrix.table.keep$loc, input_list$duplication[which(input_list$sample == patient)])
      missing.df = data.frame(duplication = missing,
                              sample = patient)
      missing.df = separate(missing.df, 
                            col = duplication,
                            into = c('Chromosome', 'Position'),
                            sep = ';',
                            remove = F)
      
      #' add artificial data for missing positions; in this case just 1
      missing.df$NOR.DP = 1
      missing.df$NOR.RD = 1
      table.out = rbind(table.out, missing.df)
    }
    locations.out = rbind(locations.out, table.out)
  }
  

  #' prepare the final output
  PON_out = as.data.frame(do.call('cbind', split(locations.out[, c('NOR.DP')], locations.out$sample)))
  row.names(PON_out) = matrix.table.keep$loc
  
  #' mean-normalization
  mean_normalization = function(x){
    x / mean(x)
  }
  
  message('Starting mean-normalization')
  
  PON_normalized = apply(PON_out, 2, mean_normalization)
  
  #' return object
  return(list(selcted_bins = matrix.table.keep$loc,
              PON_normalized = PON_normalized))
}

PON_out = prepare_PON(normal_samples = BRCA_PON_df,
                      sample_threshold = 0.95)

saveRDS(object = PON_normalized, file = 'PON_BRCA/PON_normalized.rds')

PON = readRDS('PON_BRCA/PON_normalized.rds')

#' prepare table for DryClean function:
#' we need to create a gRanges object (similar to output from fragCounter)
#' afterwards rPCA decomposition is done on matrix.

modify_PON = function(data, path_to_save){
  
  PON_path = data.table::data.table()
  normalized_data = as.data.frame(data)
  
  bins = row.names(normalized_data)
  
  for(i in 1:length(normalized_data)){
    print(i)
    sample_selected = data.frame(normalized_data[, i])
    sample_selected$seq = bins
    sample_selected_ext = separate(sample_selected, 
                                   col = seq,
                                   into = c('seqnames', 'ranges'),
                                   sep = ';')
    colnames(sample_selected_ext)[1] = 'reads.corrected'
    sample_selected_ext$start = sample_selected_ext$ranges
    sample_selected_ext$end = sample_selected_ext$ranges
    sample_selected_ext$ranges = NULL
    GR_sample = makeGRangesFromDataFrame(df = sample_selected_ext, keep.extra.columns = T)
    
    #' append one nucleotide, ot have a proper range object
    GR_sample = resize(GR_sample, width(GR_sample) + 1, fix = 'start')
    saveRDS(GR_sample, file = paste0(path_to_save, 'sample', i, '.rds'))
    
    #' prepare data table for subsequent follow-up
    paths = data.table::data.table(sample = paste0('sample', i),
                                   normal_cov = paste0('~/Documents/MSKCC/07_FacetsReview/PON_BRCA/sample', i, '.rds'))
    PON_path = rbind(PON_path, paths)
  }
  saveRDS(PON_path, file = 'PON_BRCA/normal_table.rds')
}

modify_PON(data = PON, path_to_save = 'PON_BRCA/')


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' start prepare_detergent
library(dplyr)
library(data.table)
library(pbmcapply)
library(parallel)
library(dplyr)
library(gUtils)
library(plyr)

normal_PON = readRDS('PON_BRCA/normal_table.rds')
for(i in seq(50, 500, 50)){
  table_subset = normal_PON[1:i, ]
  filename = paste0('PON_BRCA/table_', i, '.rds')
  saveRDS(object = table_subset, file = filename)
  print(filename)
  detergent = prepare_detergent(normal.table.path = filename,
                                path.to.save = 'PON_BRCA',
                                save.pon = T)
  rm(filename)
  rm(detergent)
  rm(table_subset)
}




detergent = prepare_detergent(normal.table.path = 'PON_BRCA/normal_table.rds', 
                              path.to.save = 'PON_BRCA', 
                              save.pon = T)

saveRDS(detergent, file = 'PON_BRCA/detergent_compressed.rds', compress = T)

#' out