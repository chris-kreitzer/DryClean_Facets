## Running decompostion (DryClean) on tumor samples
## start_wash_cycle()
## 
## Start: 09/18/2021
## Revision: 02/26/2022
## chris-kreitzer


clean()
library(GenomicRanges)
setwd('~/Documents/GitHub/DryClean_Facets/')


Tumor_samples = list.files('~/Documents/MSKCC/07_FacetsReview/DryClean/TUMOR_BRCA/', full.names = T)

test = Tumor_samples[1:2]

for(i in test){
  cov = readRDS(i)
  cov_converted = GRanges(as.data.frame(GRangesList(cov)))
  cov_out = dryclean::start_wash_cycle(cov = cov_converted,
                                       detergent.pon.path = '~/Documents/GitHub/DryClean_Facets/detergent.rds',
                                       whole_genome = F,
                                       chr = NA,
                                       germline.filter = FALSE,
                                       mc.cores = 4)
  saveRDS(cov_out, file = paste0('~/Desktop/', basename(i), '.rds'))
}


##-----------------------------------------------------------------------------
## Cluster Version

require('dryclean', '/juno/home/kreitzec/R/x86_64-pc-linux-gnu-library/4.0/')
require('data.table', '/juno/home/kreitzec/R/x86_64-pc-linux-gnu-library/4.0/')
require('GenomicRanges', '/juno/home/kreitzec/R/x86_64-pc-linux-gnu-library/4.0/')
library(parallel)

Tumor_samples = list.files(path = '/home/kreitzec/DryClean/v2/GrTumorNorm/TUMOR_BRCA', pattern = '^sample.*',
                           full.names = T)

for(i in Tumor_samples){ 
  cov = readRDS(i)
  cov_converted = GRanges(as.data.frame(GRangesList(cov)))
  cov_out = dryclean::start_wash_cycle(cov = cov_converted,
                                       detergent.pon.path = '/home/kreitzec/DryClean/v2/PON/detergent.rds',
                                       whole_genome = F,
                                       chr = NA,
                                       germline.filter = FALSE,
                                       mc.cores = 4)
  saveRDS(cov_out, file = paste0('/home/kreitzec/DryClean/v2/TumorClean/', basename(i), '.rds'))
}
