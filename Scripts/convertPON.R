Sys.setenv('R_MAX_VSIZE'=32000000000)
rm(list = ls())
gc()

library(tidyverse)
library(data.table)
library(GenomicRanges)
library(tidyr)
library(IRanges)




PON = read.csv('~/Documents/MSKCC/07_FacetsReview/DryClean/Data4Analysis/normalizedPON.txt', sep = '\t')



#' prepare table for DryClean function:
#' we need to create a gRanges object (similar to output from fragCounter)
#' afterwards rPCA decomposition is done on matrix.


modify_PON = function(data, path_to_save){
  tryCatch({
    message('TryCatch is running')
    PON_path = data.table::data.table()
    normalized_data = as.data.frame(data)
    row.names(normalized_data) = normalized_data$duplication
    normalized_data$duplication = NULL
    
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
                                     normal_cov = paste0('~/Documents/MSKCC/07_FacetsReview/DryClean/PON_BRCA/sample', i, '.rds'))
      PON_path = rbind(PON_path, paths)
    }
    
    saveRDS(PON_path, file = '~/Documents/MSKCC/07_FacetsReview/DryClean/PON_BRCA/normal_table.rds')
   
  },
  error = function(cond){
    message(paste('Sample: ', sample, ' failed'))
    message(cond)
    return(NA)
  })
}

modify_PON(data = PON, 
           path_to_save = '~/Documents/MSKCC/07_FacetsReview/DryClean/PON_BRCA/')



