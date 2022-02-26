## Convert Tumors to GRanges object; 
## making ready for DryCleans decompostion


Sys.setenv('R_MAX_VSIZE' = 32000000000)
rm(list = ls())
gc()

library(tidyverse)
library(data.table)
library(GenomicRanges)
library(tidyr)
library(IRanges)


Tumors = readRDS('~/Documents/GitHub/DryClean_Facets/Data_out/TumorNormalizedAll.rds')



#' prepare table for DryClean decompostion:
#' we need to create a gRanges object (similar to output from fragCounter)
#' afterwards rPCA decomposition is done on matrix.


modify_Tumors = function(data, 
                         path_to_save){
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
                         tumor_cov = paste0(path_to_save, 'sample', i, '.rds'))
      PON_path = rbind(PON_path, paths)
    }
    
    saveRDS(PON_path, file = paste0(path_to_save, 'tumor_table.rds'))
    
  },
  error = function(cond){
    message(paste('Sample: ', sample, ' failed'))
    message(cond)
    return(NA)
  })
}

modify_Tumors(data = Tumors,
              path_to_save = '~/Documents/MSKCC/07_FacetsReview/DryClean/TUMOR_BRCA/')

#' out
y = readRDS('~/Documents/MSKCC/07_FacetsReview/DryClean/TUMOR_BRCA/tumor_table.rds')
a = readRDS('~/Documents/MSKCC/07_FacetsReview/DryClean/TUMOR_BRCA/sample1.rds')
a = annoGR2DF(a)
View(a)
