## Identify germline events from the PON
## For this purpose, normal samples are treated like tumors:
## 
## Firstly, we need to decompse the normals, then we will run the
## identify_gmerline() function
## importantly - normals are treated like tumors

## start: 09/17/2021
## revision: 02/22/2022
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
normal_samples_selected = sample(x = normal_samples, size = 15, replace = T)
dir.create(path = paste0(PON_path, 'decomposed_samples'))

for(i in 1:length(normal_samples_selected)){
  decomp.1 = start_wash_cycle(cov = normal_samples_selected[i], 
                              detergent.pon.path = "detergent.rds", 
                              whole_genome = FALSE, 
                              chr = NA, 
                              germline.filter = FALSE, 
                              mc.cores = 6)
  saveRDS(object = decomp.1, file = paste0(PON_path, 'decomposed_samples/', basename(normal_samples_selected[i])))
  rm(decomp.1)
  gc()
}







decomp.1 = start_wash_cycle(cov = sample.1, detergent.pon.path = "~/git/dryclean/inst/extdata/", whole_genome = TRUE, chr = NA, germline.filter = FALSE)



working_path = ' '

Germline = identify_germline(normal.table.path = paste0(working_path, 'decomposed_samples/normal_table.rds'), 
                             path.to.save = paste0(working_path, 'decomposed_samples'), 
                             signal.thresh = 0.5, 
                             pct.thresh = 0.98)

decomp.1 = start_wash_cycle(cov = '~/Documents/MSKCC/07_FacetsReview/DryClean/PON_BRCA/sample1.rds', 
                            detergent.pon.path = "detergent.rds", 
                            whole_genome = T, chr = NA, germline.filter = FALSE)

germ1 = identify_germline()

j = readRDS('~/Documents/MSKCC/07_FacetsReview/DryClean/PON_BRCA/sample1.rds')
j
