## The third part in the DryClean process;
## modify the tumor samples and process like described
## 
## 09/11/2021
## chris-kreitzer



#' modify the tumor samples; work on total read depth first;
#' we need a coverage file alike we needed for the normal


## Libraries, Dependencies and Input
# library(tidyverse)
# library(vroom)
# library(GenomicRanges)
Sys.setenv("VROOM_SHOW_PROGRESS" = "false")
require('tidyverse', '/juno/home/kreitzec/R/x86_64-pc-linux-gnu-library/4.0/')
require('vroom', '/juno/home/kreitzec/R/x86_64-pc-linux-gnu-library/4.0/')
require('GenomicRanges', '/juno/home/kreitzec/R/x86_64-pc-linux-gnu-library/4.0/')


## INPUT:
#' sample_path = absolute path to dataframe which contains the paths for respective tumor samples
#' PON_path = one example path which contains the positions that were used for the PON creation (.rds format)
#' path_to_save = where should the converted tumor files be stored
#' threshold = minimum position overlap, between PON positions and tumor samples to be considered


#' prepare tumor:
prepare_tumor_array = function(sample_path, 
                               PON_path, 
                               path_to_save, 
                               threshold = NULL){
  
  message('Please provide a path to a normal (PON) sample [.rds]')
  message('Tumor list must contain "sample -column" which is path to count__matrix')
  
  if(is.null(threshold)){
    stop('process is interrupted with no position-threshold', call. = T)
  }
  
  #' marker positions which are used for the creation of the PON
  PON_used = readRDS(PON_path)
  Markers_used = as.data.frame(PON_used)
  Markers_used = data.frame(chromosome = Markers_used$seqnames,
                            position = Markers_used$start)
  Markers_used$merged_position = paste(Markers_used$chromosome, Markers_used$position, sep = ';')
  
  #' import tumor samples list which is to be analyzed
  tumor_list = read.csv(sample_path, sep = '\t')
  tumor_list = data.frame(sample = tumor_list[1:2, ])
  
  extracted_tumors = data.frame()
  NA_tumors = c()
  for(i in 1:nrow(tumor_list)){
    data.in = vroom::vroom(tumor_list$sample[i])
    data.in = data.frame(Chromosome = data.in$Chromosome,
                         Position = data.in$Position,
                         depth = data.in$File2R + data.in$File2A)
    data.in$sample = basename(tumor_list$sample[i])
    data.in$duplication = paste(data.in$Chromosome, data.in$Position, sep = ';')
    
    try({
      tumor_sample_out = data.in[which(data.in$duplication %in% Markers_used$merged_position), ]
      if(dim(tumor_sample_out)[1] / dim(Markers_used)[1] <= threshold){
        print(paste0('Sample ', i, ' cannot be used for estimation, because too little bins'))
        NA_tumors = c(NA_tumors, i)
        next
        
      } else {
        missing = setdiff(Markers_used$merged_position, tumor_sample_out$duplication)
        
        if(length(missing) == 0){
          extracted_tumors = rbind(extracted_tumors, tumor_sample_out)
          
        } else {
          missing.df = data.frame(duplication = missing,
                                  sample = basename(tumor_list$sample[i]))
          missing.df = separate(missing.df,
                                col = duplication,
                                into = c('Chromosome', 'Position'),
                                sep = ';',
                                remove = F)
          missing.df$depth = 1
          extracted_tumors = rbind(extracted_tumors, tumor_sample_out, missing.df)
          
        }
      }
    })
  } 
  
  saveRDS(extracted_tumors, file = paste0(path_to_save, 'tumors_all.rds'))
  
  #' make the mean normalization, create GRobject and save the output for analysis
  for(i in unique(extracted_tumors$sample)){
    data = extracted_tumors[which(extracted_tumors$sample == i), ]
    data$normalized.depth = data$depth / mean(data$depth)
    data = data[, c('Chromosome', 'Position', 'normalized.depth', 'sample')]
    Tumor_GR = data.frame(seqnames = data$Chromosome,
                          start = as.numeric(data$Position),
                          end = as.numeric(data$Position) + 1,
                          reads.corrected = data$normalized.depth)
    Tumor_GRobject = makeGRangesFromDataFrame(df = Tumor_GR, keep.extra.columns = T)
    saveRDS(Tumor_GRobject, file = paste0(path_to_save, i, '.rds'))
    
  }
}
  
prepare_tumor_array(sample_path = '~/Desktop/test.txt', PON_path = '~/Documents/MSKCC/07_FacetsReview/PON_BRCA/sample1.rds', path_to_save = '~/Desktop/', threshold = 0.95)




## TEST run
#' INPUT:
# sample_path = data.frame(sample = c('~/Desktop/mnt/ATMcountdata/countsMerged____P-0002273-T01-IM3_P-0002273-N01-IM3.dat.gz',
#                 '~/Desktop/mnt/ATMcountdata/countsMerged____P-0003139-T02-IM5_P-0003139-N01-IM5.dat.gz'))
# write.table(sample_path, file = '~/Desktop/Test_tumor.txt', sep = '\t', row.names = F)
# 
# prepare_tumor_array(sample_path = '~/Desktop/Test_tumor.txt', PON_path = '~/Documents/MSKCC/07_FacetsReview/PON_BRCA/sample1.rds', 
#                     path_to_save = '~/Desktop/', threshold = 0.65)


#' out
  