## The third part in the DryClean process;
## modify the tumor samples and process like described
## 
## 
## Germline is running
## 



#' modify the tumor samples; work on total read depth first;
#' we need a coverage file alike we needed for the normal


## Libraries, Dependencies and Input
Sys.setenv("VROOM_SHOW_PROGRESS" = "false")
library(tidyverse)


#' prepare tumor:
prepare_tumor_array = function(sample_path, PON_path, threshold = NULL){
  message('Please provide a path to a normal (PON) sample [.rds]')
  message('tumor list should contain the column SAMPLE which has the path to sample')
  
  #' marker positions which are used for the creation of the PON
  PON_used = readRDS(PON_path)
  Markers_used = as.data.frame(PON_used)
  Markers_used = data.frame(chromosome = Markers_used$seqnames,
                            position = Markers_used$start)
  Markers_used$merged_position = paste(Markers_used$chromosome, Markers_used$position, sep = ';')
  
  #' import tumor samples list which is to be analyzed
  tumor_list = read.csv(sample_path, sep = '\t')
  
  #' think about lapply alternative
  data_merged = data.frame()
  for(i in 1:nrow(tumor_list)){
    data.in = vroom::vroom(tumor_list$sample[i], show_col_types = FALSE)
    data.in = data.frame(Chromosome = data.in$Chromosome,
                         Position = data.in$Position,
                         depth = data.in$File2R + data.in$File2A)
    data.in$sample = basename(tumor_list$sample[i])
    data.in = data.in[which(data.in$depth > 30), ]
    
    # stopifnot(ncol(data.in) != 9)
    data_merged = rbind(data_merged, data.in)
  }
  
  #' check for duplicated entries (chromosome; position)
  data_merged$duplication = paste(data_merged$Chromosome, data_merged$Position, sep = ';')
  container = c()
  for(i in unique(data_merged$sample)){
    if(any(duplicated(data_merged$duplication[which(data_merged$sample == i)]))){
      container = c(container, i)
    } 
    else next
  }
  
  message(paste0(length(container), ' samples have duplicated bins'), appendLF = T)
  data_merged = data_merged[!data_merged$sample %in% container,, drop = F]
  rm(container)
  message('Quality control ended')
  message('Welcome. I am creating an n x m matrix which is equal among all normal samples')
  
  #' extract respective positions from tumor samples and create GRange object
  extracted_tumors = data.frame()
  NA_tumors = c()
  for(i in unique(data_merged$sample)){
    tumor_sample_out = data_merged[which(data_merged$sample == i & data_merged$duplication %in% Markers_used$merged_position), ]
    if(dim(tumor_sample_out)[1] / dim(Markers_used)[1] <= 0.95){
      print(paste0('Sample ', i, 'cannot be used for estimation, because too little bins'))
      NA_tumors = c(NA_tumors, i)
      next
    } else if (dim(tumor_sample_out)[1] / dim(Markers_used)[1] > 0.95){
      missing = setdiff(Markers_used$merged_position, data_merged$duplication[which(data_merged$sample == i)])
      if(length(missing) == 0){
        extracted_tumors = rbind(extracted_tumors, tumor_sample_out)
      } else {
        missing.df = data.frame(duplication = missing,
                                sample = i)
        missing.df = separate(missing.df,
                              col = duplication,
                              into = c('Chromosome', 'Position'),
                              sep = ';',
                              remove = F)
        missing.df$depth = 1
        extracted_tumors = rbind(extracted_tumors, tumor_sample_out)
      }
    } 
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



coverage_file = readRDS("~/git/dryclean/data/dummy_coverage.rds")
coverage_file