## Running dryclean on tumor sample within R
## Firstly, we need to modify the input (countmatrix to GRange object)
##
## 
## start: 09/11/2021
## revision: 02/23/2022
## chris-kreitzer


clean()
setup(working.path = '~/Documents/GitHub/DryClean_Facets/')




#' Modify the tumor samples; 
#' We need the same bins (dimensions) for the tumors alike the PON
prepareTumors = function(sample_path,
                         PON_path,
                         path_to_save){
  
  #' reference bins (from PON)
  PON = readRDS(PON_path)
  bin_range = paste(seqnames(PON), ranges(PON), sep = ';')
  ref_bin = gsub(pattern = '\\-.*', replacement = '', x = bin_range)
  
  
  #' structure checks
  stopifnot(is(sample_path, 'data.table'))
  stopifnot(is(PON, 'GRanges'))  
  message('Only positions >= 35 reads are supported')
  
  FLAG = c()
  
  #' remove samples with duplicated entries (faulty snp-pileup)
  mat_tumor = mclapply(sample_path$sample, function(x_tumor){
    print(x_tumor)
    tumor = tryCatch(data.table::fread(x_tumor), error = function(e) NULL)
    sample_name = substr(x = basename(x_tumor), start = 17, stop = 33)
  
    if(!is.null(tumor)){
      tumor$depth = tumor$File2A + tumor$File2R
      tumor$sample = sample_name
      tumor = tumor[which(tumor$depth >= 35), c('Chromosome', 'Position', 'depth')]
      tumor$duplication = paste(tumor$Chromosome, tumor$Position, sep = ';')
      tumor$sample = sample_name
      
      if(any(duplicated(tumor$duplication))){
        FLAG = c(FLAG, sample_name)
        print(FLAG)
      }
      
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
    }
    return(Tumor_out)
  }, mc.cores = 1)
 
  gc()
  
  mat.all = rbindlist(mat_tumor, fill = T)
  
  return(mat.all)
}
  
  
  
  
test = data.table(sample = c('~/Desktop/mnt/ATMcountdata/countsMerged____P-0002273-T01-IM3_P-0002273-N01-IM3.dat.gz', 
                             '~/Desktop/mnt/ATMcountdata/countsMerged____P-0003139-T02-IM5_P-0003139-N01-IM5.dat.gz'))

x = prepareTumors(sample_path = test, PON_path = '~/Documents/MSKCC/07_FacetsReview/DryClean/PON_BRCA/sample1.rds')




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
                               save_compendium = NULL,
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
  
  extracted_tumors = data.frame()
  NA_tumors = c()
  for(i in 1:nrow(tumor_list)){
    data.in = vroom::vroom(tumor_list$sample[i])
    data.in = data.frame(Chromosome = data.in$Chromosome,
                         Position = data.in$Position,
                         depth = data.in$File2R + data.in$File2A)
    data.in = data.in[which(data.in$depth > 30), ]
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
  
  #' save a comprehensive file if required
  if(!is.null(save_compendium)){
    saveRDS(extracted_tumors, file = paste0(path_to_save, 'tumors_all.rds'))  
  }
  
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
  
prepare_tumor_array(sample_path = '~/Desktop/test.txt', 
                    PON_path = '~/Documents/MSKCC/07_FacetsReview/PON_BRCA/sample1.rds', 
                    path_to_save = '~/Desktop/', 
                    threshold = 0.95)



## TEST run
#' INPUT:
# sample_path = data.frame(sample = c('~/Desktop/mnt/ATMcountdata/countsMerged____P-0002273-T01-IM3_P-0002273-N01-IM3.dat.gz',
#                 '~/Desktop/mnt/ATMcountdata/countsMerged____P-0003139-T02-IM5_P-0003139-N01-IM5.dat.gz'))
# write.table(sample_path, file = '~/Desktop/Test_tumor.txt', sep = '\t', row.names = F)
# 
# prepare_tumor_array(sample_path = '~/Desktop/Test_tumor.txt', PON_path = '~/Documents/MSKCC/07_FacetsReview/PON_BRCA/sample1.rds', 
#                     path_to_save = '~/Desktop/', threshold = 0.65)


#' out





## Quality Control of the conversion;
Tumor_sample_in = vroom::vroom('countsMerged____P-0000584-T03-IM6_P-0000584-N01-IM6.dat.gz')
PON = readRDS('~/Documents/MSKCC/07_FacetsReview/PON_BRCA/sample1.rds')

a = Tumor_sample_in[which(Tumor_sample_in$Position == '1001150' & Tumor_sample_in$Chromosome == 1), ]
b = Tumor_sample_in[which(Tumor_sample_in$Position == '2487984' & Tumor_sample_in$Chromosome == 1), ]
c = Tumor_sample_in[which(Tumor_sample_in$Position == '152861750' & Tumor_sample_in$Chromosome == 'X'), ]
d = Tumor_sample_in[which(Tumor_sample_in$Position == '29910760' & Tumor_sample_in$Chromosome == 6), ]

#'
Tumor_sample = data.frame(Chromosome = Tumor_sample_in$Chromosome,
                     Position = Tumor_sample_in$Position,
                     depth = Tumor_sample_in$File2R + Tumor_sample_in$File2A)
Tumor_sample = Tumor_sample[which(Tumor_sample$depth > 30), ]
Tumor_sample$duplication = paste(Tumor_sample$Chromosome, Tumor_sample$Position, sep = ';')

#'
PON_used = PON
Markers_used = as.data.frame(PON_used)
Markers_used = data.frame(chromosome = Markers_used$seqnames,
                          position = Markers_used$start)
Markers_used$merged_position = paste(Markers_used$chromosome, Markers_used$position, sep = ';')

#'
out = Tumor_sample[which(Tumor_sample$duplication %in% Markers_used$merged_position), ]

#'
missing = setdiff(Markers_used$merged_position, out$duplication)
missing.df = data.frame(duplication = missing)
missing.df = separate(missing.df,
                      col = duplication,
                      into = c('Chromosome', 'Position'),
                      sep = ';',
                      remove = F)
missing.df$depth = 1
out_full = rbind(out, missing.df)

#' 
out_full$normalized = out_full$depth / mean(out_full$depth)
plot(density(out_full$normalized))
#' the mean is centered around 1

#' out





  