##############################
## Tumor_countmatrix_conversion
##############################
#' @name prepareTumors
#'
#' @description: Adjust tumor count-matrices (snp-pileup), to 
#' appropriate dimension (equally like PON). This function substitute
#' missing bins (PON reference), and conducts the mean normalization
#' 
#' @export
#' @param sample_path data.table(); path to tumor-matrices (absolute path); column sample = absolute path
#' @param PON_path rds(); path to reference PON (required to extract respective bins); .rds object
#' @param path_to_save path where the return dataframe() will be saved
#' 
#' @return dataframe()
#' @author chris-kreitzer

## start: 09/11/2021
## revision: 02/23/2022
## revision: 04/13/2022


clean()
setup(working.path = '~/Documents/GitHub/DryClean_Facets/')

library(GenomicRanges)
library(IRanges)
library(GenomeInfoDb)
library(parallel)
library(data.table)

#' Adjust tumor count-matrices
#' Equal dimensions (bins) like PON
prepareTumors = function(sample_path,
                         PON_path,
                         path_to_save,
                         mc.cores = 1){
  
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
      Tumor_out = Tumor_out[, c(4,6)]
    }
    return(Tumor_out)
  }, mc.cores = mc.cores)
 
  gc()
  
  mat.all = Reduce(function(...) merge(..., all = TRUE, by = 'duplication'), mat_tumor)
  
  return(mat.all)
  
  message('The output will not be saved automatically')
  
}
  

#' example: 
test = data.table(sample = c('~/Desktop/mnt/ATMcountdata/countsMerged____P-0002273-T01-IM3_P-0002273-N01-IM3.dat.gz', 
                             '~/Desktop/mnt/ATMcountdata/countsMerged____P-0003139-T02-IM5_P-0003139-N01-IM5.dat.gz'))

x = prepareTumors(sample_path = test, PON_path = '~/Documents/MSKCC/07_FacetsReview/DryClean/PON_BRCA/sample1.rds')





##-----------------------------------------------------------------------------
## Cluster Version
## Libraries, Dependencies and Input
# library(tidyverse)
# library(vroom)
# library(GenomicRanges)
require('tidyverse', '/juno/home/kreitzec/R/x86_64-pc-linux-gnu-library/4.0/')
require('data.table', '/juno/home/kreitzec/R/x86_64-pc-linux-gnu-library/4.0/')
require('GenomicRanges', '/juno/home/kreitzec/R/x86_64-pc-linux-gnu-library/4.0/')
library(parallel)

## INPUT:
#' sample_path = data.table() with ABSOLUTE path to tumor count_matrices
#' PON_path = one example path which contains the positions that were used for the PON creation (.rds format)
#' path_to_save = where should the converted tumor files be stored
Tumor_samples = list.files(path = '/juno/home/kreitzec/DryClean/v2/Tumor_CountFiles/', full.names = T)
Tumor_samples = data.table(sample = Tumor_samples)
message(paste0('There are: ', nrow(Tumor_samples), ' available for analysis'))

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
    tumor = tryCatch(data.table::fread(x_tumor), error = function(e) NULL)
    sample_name = substr(x = basename(x_tumor), start = 17, stop = 33)
    print(sample_name)
    
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
      Tumor_out = Tumor_out[, c(4,6)]
    }
    return(Tumor_out)
  }, mc.cores = 1)
  
  gc()
  
  message('Objects are merged now')
  mat.all = Reduce(function(...) merge(..., all = TRUE, by = 'duplication'), mat_tumor)
  
  saveRDS(object = mat.all, file = '/juno/home/kreitzec/DryClean/v2/TumorNormalizedAll.rds')
  
}

prepareTumors(sample_path = Tumor_samples, 
              PON_path = '/home/kreitzec/DryClean/v2/PON/PON_BRCA/sample1.rds')




