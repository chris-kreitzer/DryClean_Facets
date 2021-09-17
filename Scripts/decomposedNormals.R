## The second step in DryClean is to identify germline events,
## in order to remove them from the panel of normals;
##
## 09/06/2021
## chris-kreitzer


rm(list = ls())
.rs.restartR()
setwd('~/Documents/GitHub/DryClean_Facets/')
set.seed(12)


## Libraries and Input
library(gUtils)
library(dryclean)


## Processing
working_path = c('~/Documents/MSKCC/07_FacetsReview/PON_BRCA/')
PON_table = readRDS(paste0(working_path, 'normal_table.rds'))
PON_table$decomposed_cov = NA
dir.create(path = paste0(working_path, 'decomposed_samples'))

message(paste0('The decomposition of the germline consits of n = ', nrow(PON_table), ' samples'))

for(i in 1:nrow(PON_table)){
  print(i)
  decomposed_sample = dryclean::start_wash_cycle(cov = PON_table$normal_cov[i], 
                                                 detergent.pon.path = paste0(working_path, 'detergent.rds'), 
                                                 whole_genome = F)
  saveRDS(decomposed_sample, file = paste0(working_path, 'decomposed_samples/', basename(PON_table$normal_cov[i])))
  PON_table$decomposed_cov[i] = paste0(working_path, 'decomposed_samples/', basename(PON_table$normal_cov[i]))
  saveRDS(PON_table, file = paste0(working_path, 'decomposed_samples/normal_table.rds'))
}


#' out
