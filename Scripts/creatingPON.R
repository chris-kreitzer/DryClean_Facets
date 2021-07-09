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
    print(chromosome)
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
  
  #' prepare the final output
  PON_out = as.data.frame(do.call('cbind', split(sample_PON[, c('NOR.DP')], sample_PON$sample)))
  row.names(PON_out) = paste(bins_PON$chromosome, bins_PON$loc, sep = ';')
  
  #' mean-normalization
  mean_normalization = function(x){
    (x - mean(x)) / (max(x) - min(x))
  }
  
  message('Starting mean-normalization')
  
  PON_normalized = apply(PON_out, 2, mean_normalization)
  
  #' return object
  return(list(selcted_bins = bins_PON,
              PON_out = PON_out,
              PON_normalized = PON_normalized))
}

PON_out = prepare_PON(normal_samples = BRCA_PON_df)



#' prepare table for DryClean function:
#' we need to create a gRanges object (similar to output from fragCounter)
#' afterwards rPCA decomposition is done on matrix.

modify_PON = function(data, is_list_input = NULL, path_to_save){
  
  if(!is.null(is_list_input)){
    normalized_data = as.data.frame(data)
  } else {
    normalized_data = as.data.frame(data$PON_normalized)
  }
  
  bins = row.names(normalized_data)
  normalized_data = as.data.frame(normalized_data)
  
  for(i in 1:length(normalized_data)){
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
    saveRDS(GR_sample, file = paste0(path_to_save, 'sample', i, '.rds'))
  }
}

u = PON_out$PON_normalized
u = u[, 1:5]

modify_PON(data = u, is_list_input = F, path_to_save = 'PON_BRCA/')



head(u)
library(tidyverse)

#' convert dataframe into GRanges object

#' example 1
x = PON_out$PON_normalized
row_names = row.names(x)
x = as.data.frame(x)
x = data.frame(x$`P-0004835-T02-IM6`)
x$seq = row_names

x = separate(data = x, col = seq, into = c('seqnames', 'ranges'), sep = ';')
colnames(x)[1] = 'reads.corrected'
x$start = x$ranges
x$end = x$ranges
x$ranges = NULL

y = makeGRangesFromDataFrame(df = x, keep.extra.columns = T)

saveRDS(y, file = 'PON_BRCA/sample1.rds')

#' example2
x = PON_out$PON_normalized
row_names = row.names(x)
x = as.data.frame(x)
x = data.frame(x$`P-0014407-T01-IM6`)
x$seq = row_names
x = separate(data = x, col = seq, into = c('seqnames', 'ranges'), sep = ';')
colnames(x)[1] = 'reads.corrected'
x$start = x$ranges
x$end = x$ranges
x$ranges = NULL

y = makeGRangesFromDataFrame(df = x, keep.extra.columns = T)
saveRDS(y, file = 'PON_BRCA/sample2.rds')


#' prepare the table:
normal_path = data.table::data.table(sample = c('sample1', 'sample2'),
                         normal_cov = c('~/Documents/MSKCC/07_FacetsReview/PON_BRCA/sample1.rds',
                                        '~/Documents/MSKCC/07_FacetsReview/PON_BRCA/sample2.rds'))


saveRDS(normal_path, file = 'PON_BRCA/normal_table.rds')

detergent = dryclean::prepare_detergent(normal.table.path = 'PON_BRCA/normal_table.rds', 
                                        use.all = T, 
                                        choose.randomly = F, path.to.save = 'PON_BRCA', 
                                        save.pon = T)



a = readRDS('PON_BRCA/detergent.rds')
a$U.hat



#' make a cross-check if matrix is okay;

#' start prepare_detergent - DryClean (and documentation)









