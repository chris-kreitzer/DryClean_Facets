## Identify germline events from the PON
## For this purpose, normal samples are treated like tumors:
## 
## Firstly, we need to decompse the normals, then we will run the
## identify_gmerline() function
## importantly - normals are treated like tumors

## start: 09/17/2021
## revision: 02/22/2022
## revision out: Addy suggested that this function should be taken down!
## chris-kreitzer

clean()
setup()

library(parallel)
library(dryclean)


## 2. Identifying germline events
#' Since a PON is used for decomposing the tumor samples, a method is required 
#' identify and remove germline events. This is achieved by looking at all normal 
#' samples as a population and infer the markers that have a copynumber events at 
#' a given frequency, set by user.

## Decomposition of the Normals with the start_wash_cycle() function
## I will run this function on 15 randomly selected normals (due to time issues)
PON_path = '~/Documents/MSKCC/07_FacetsReview/DryClean/PON_BRCA/'
normal_samples = list.files(path = PON_path, pattern = '^sample.*', full.names = T)
normal_samples_selected_x = sample(x = normal_samples, size = 15, replace = T)
dir.create(path = paste0(PON_path, 'decomposed_samples'))

for(i in 1:length(normal_samples_selected_x)){
  decomp.1 = start_wash_cycle(cov = normal_samples_selected_x[i], 
                              detergent.pon.path = "detergent.rds", 
                              whole_genome = FALSE, 
                              chr = NA, 
                              germline.filter = FALSE, 
                              mc.cores = 6)
  saveRDS(object = decomp.1, file = paste0(PON_path, 'decomposed_samples/', basename(normal_samples_selected_x[i])))
  rm(decomp.1)
  gc()
}


#' prepare full table for germline identification:
decomposed_samples = list.files(path = '~/Documents/MSKCC/07_FacetsReview/DryClean/PON_BRCA/decomposed_samples/', full.names = T)

normal_table = data.frame(decomposed_cov = decomposed_samples)
normal_table$normal_cov = paste0('/Users/chriskreitzer/Documents/MSKCC/07_FacetsReview/DryClean/PON_BRCA/', basename(decomposed_samples))
normal_table$sample = basename(normal_table$decomposed_cov)
normal_table = normal_table[, c(3,2,1)]
normal_table = as.data.table(normal_table)
saveRDS(normal_table, file = '~/Documents/MSKCC/07_FacetsReview/DryClean/PON_BRCA/normal_table_n30.rds')


## identify Germline events;
Germline = identify_germline(normal.table.path = '~/Documents/MSKCC/07_FacetsReview/DryClean/PON_BRCA/normal_table_n30.rds', 
                             path.to.save = '~/Documents/MSKCC/07_FacetsReview/DryClean/PON_BRCA/decomposed_samples/', 
                             signal.thresh = 0.5, 
                             pct.thresh = 0.98)

