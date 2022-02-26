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
      
      #' append one nc, to have a proper range object
      GR_sample = resize(GR_sample, width(GR_sample) + 1, fix = 'start')
      saveRDS(GR_sample, file = paste0(path_to_save, 'sample', i, '.rds'))
      
      #' prepare data table for subsequent follow-up
      paths = data.table::data.table(sample = paste0('sample', i),
                                     tumor_cov = paste0(path_to_save, 'sample', i, '.rds'))
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

modify_Tumors(data = Tumors,
              path_to_save = '~/Documents/MSKCC/07_FacetsReview/DryClean/TUMOR_BRCA/')

#' out








clean()
gc()


x = readRDS('~/Documents/GitHub/DryClean_Facets/Data_out/TumorNormalizedAll.rds')
PON = readRDS('~/Documents/MSKCC/07_FacetsReview/DryClean/PON_BRCA/sample1.rds')
ref_bins = paste(seqnames(PON), ranges(PON), sep = ';')
ref_bin = gsub(pattern = '\\-.*', replacement = '', x = ref_bins)



sample1 = data.table::fread('~/Desktop/countsMerged____P-0002713-T02-IM6_P-0002713-N01-IM6.dat.gz')
sample1$depth = sample1$File2A + sample1$File2R
sample1 = sample1[which(sample1$depth >= 35), ]
sample1 = sample1[,c(1,2,13)]
sample1$duplication = paste(sample1$Chromosome, sample1$Position, sep = ';')


tumor = sample1
sample_name = 'P-0002713'
tumor$sample = sample_name

#' substitute missing ref_bins
missing_bins = setdiff(ref_bin, tumor$duplication)
missing_df = data.table(Chromosome = unlist(strsplit(as.character(missing_bins), ';'))[2*(1:length(missing_bins)) - 1],
                        Position = unlist(strsplit(as.character(missing_bins), ';'))[2*(1:length(missing_bins))],
                        duplication = missing_bins, 
                        sample = sample_name,
                        depth = 1)

Tumor_out = rbind(tumor[which(tumor$duplication %in% ref_bin), ],
                  missing_df)

#' mean normalization
Tumor_out$norm_mean = NA
Tumor_out$norm_mean[which(Tumor_out$depth != 1)] = Tumor_out$depth[which(Tumor_out$depth != 1)] / mean(Tumor_out$depth[which(Tumor_out$depth != 1)])
norm.mean = mean(Tumor_out$norm_mean[which(Tumor_out$depth != 1)])
Tumor_out$norm_mean[which(Tumor_out$depth == 1)] = norm.mean
colnames(Tumor_out)[ncol(Tumor_out)] = paste(sample_name, colnames(Tumor_out)[ncol(Tumor_out)], sep = ';')
Tumor_out = Tumor_out[, c(4,6)]
}
return(Tumor_out)
