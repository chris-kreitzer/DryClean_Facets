##############################
## Tumors_2_Granges
##############################
#' @name Tumors_2_Granges
#'
#' @description: This function takes the mean-normalized tumor count-
#' matrix (equal dimension then PON) and creates individual GRanges objects.
#' This step is required for the decomposition with the detergent.
#' 
#' @export
#' @param data dataframe(); comprehensive, mean-normalized tumor count-matrices
#' @param path_to_save absolute_path(); where the individual converted tumor samples should be saved
#' 
#' @return NULL. Individual tumor samples are saved in path provided above
#' @author chris-kreitzer

## start: 10/05/2021
## revision: 02/28/2022
## revision: 04/13/2022



Sys.setenv('R_MAX_VSIZE' = 32000000000)
clean()
gc()

library(tidyverse)
library(data.table)
library(GenomicRanges)
library(tidyr)
library(IRanges)

#' load the input data (normalized tumor counts)
Tumors = readRDS('~/Documents/GitHub/DryClean_Facets/Data_out/TumorNormalizedAll.rds')


Tumors_2_Granges = function(data, 
                            path_to_save){
  
  tryCatch({
    PON_path = data.frame()
    normalized_data = as.data.frame(data)
    Y_markers = grep(pattern = 'Y.*', x = normalized_data$duplication)
    normalized_data = normalized_data[-Y_markers, ]
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
      GR_sample = BiocGenerics::sort(GR_sample)
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

#' example
Tumors_2_Granges(data = Tumors,
                 path_to_save = '~/Documents/MSKCC/07_FacetsReview/DryClean/TUMOR_BRCA/')


#' out
