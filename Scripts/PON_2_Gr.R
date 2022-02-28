## Take the comprehensive PON normalized (including: bin-filling + mean_normalization)
## extract every normal sample and create GRange object (necessary for detergent preparation)
## 
## start: 09/28/2021
## revision: 02/28/2022
## chris-kreitzer



Sys.setenv('R_MAX_VSIZE' = 32000000000)
clean()
setup(working.path = '~/Documents/GitHub/DryClean_Facets/')
gc()


library(tidyverse)
library(data.table)
library(GenomicRanges)
library(tidyr)
library(IRanges)

#' Load input dataframe
PON = data.table::fread('~/Documents/MSKCC/07_FacetsReview/DryClean/Data4Analysis/normalizedPON.txt', sep = '\t')


#' we need to create a GRanges objects for every individual NORMAL
#' afterwards rPCA decomposition is done on matrix.

modify_PON = function(data, path_to_save){
  tryCatch({
    PON_path = data.frame()
    normalized_data = as.data.frame(data)
    bins = normalized_data$duplication
    normalized_data$duplication = NULL
    
    for(i in 1:length(normalized_data)){
      sample_name = sub(pattern = ';.*$','', colnames(normalized_data)[i])
      print(sample_name)
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
      
      #' append one nc, to have a proper range object
      GR_sample = resize(GR_sample, width(GR_sample) + 1, fix = 'start')
      saveRDS(GR_sample, file = paste0(path_to_save, 'sample', i, '.rds'))
      
      #' prepare data table for subsequent follow-up
      paths = data.frame(original = sample_name,
                         sample = paste0('sample', i),
                         normal_cov = paste0(path_to_save, 'sample', i, '.rds'))
      PON_path = rbind(PON_path, paths)
    }
    
    saveRDS(PON_path, file = paste0(path_to_save, 'normal_table.rds'))
   
  },
  error = function(cond){
    message(paste('Sample: ', sample, ' failed'))
    message(cond)
    return(NA)
  })
}

modify_PON(data = PON,
           path_to_save = '~/Documents/MSKCC/07_FacetsReview/DryClean/PON_BRCA/')

#' out
