## The third part in the DryClean process;
## modify the tumor samples and process like described
## 
## 
## Germline is running
## 



#' modify the tumor samples; work on total read depth first;
#' we need a coverage file alike we needed for the normal

#' prepare tumor:
prepare_tumor_array = function(sample_path, PON_path, threshold = NULL){
  message('Please provide a path to a normal (PON) sample [.rds]')
  message('tumor list should contain the column SAMPLE which has the path to sample')
  
  #' marker positions which are used for the creation of the PON
  PON_used = readRDS(PON_path)
  Markers_used = as.data.frame(PON_used)
  Markers_used = data.frame(chromosome = Markers_used$seqnames,
                            position = Markers_used$start)
  
  #' import tumor samples list which is to be analyzed
  tumor_list = read.csv(sample_path, sep = '\t')
  
  #' create 
  for(i in unique(tumor_list$sample))
  
  lapply(y$sample, function(x) rbind.data.frame(x))
  
  
  
  
  stopifnot(ncol())
  a = readRDS('~/Documents/MSKCC/07_FacetsReview/PON_BRCA/sample1.rds')
  
  a = vroom::vroom('~/Desktop/mnt/ATMcountdata/countsMerged____P-0002273-T01-IM3_P-0002273-N01-IM3.dat.gz')
  
  x = list.files('~/Desktop/mnt/ATMcountdata/', full.names = T)
  x = x[1:2]
  y = data.frame(sample = x)
  
  
  
  
  head(a)
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



coverage_file = readRDS("~/git/dryclean/data/dummy_coverage.rds")
coverage_file