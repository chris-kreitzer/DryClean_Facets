## Running the decompostion on tumor samples
## 


library(GenomicRanges)
setwd('~/Documents/GitHub/DryClean_Facets/Tumor_samples/')
Tumor_samples = list.files(getwd(), full.names = T)

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
  saveRDS(cov_out, file = paste0('~/Desktop/cleaned_tumors/', basename(i), '.rds'))
}

