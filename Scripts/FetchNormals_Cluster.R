##-----------------------------------------------------------------------------
## The following lines are hard-coded for Cluster (xjuno) use
#' fetch coordinates from Normal samples; 
require('facets', '/juno/home/kreitzec/R/x86_64-pc-linux-gnu-library/4.0/')
PON = read.csv('/juno/home/kreitzec/DryClean/PON_BRCA_Paths.txt', sep = '\t', header = F, row.names = F)

BRCA_PON_list = list()
for(i in unique(PON)){
  data.in = facets::readSnpMatrix(i)
  data.processed = data.in[which(data.in$NOR.DP >= 35), c('Chromosome', 'Position', 'NOR.DP', 'NOR.RD')]
  data.processed$sample = substr(x = basename(i), start = 17, stop = 33)
  BRCA_PON_list[[i]] = data.processed
}

saveRDS(object = BRCA_PON_list, file = '/juno/home/kreitze/DryClean/BRCA_PON_list.rds')

#' out