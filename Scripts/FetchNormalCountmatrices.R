## This function will select and merge x NORMAL samples
##
## This script is intended to run on the xjuno cluster
## BRCA data (gene panel 468) was fetched from cBIO: 07/07/2021
## 
## start: 07/18/2021
## revision: 02/28/2022
## chris-kreitzer

clean()
setup(working.path = '~/Documents/GitHub/DryClean_Facets/')

## Input data
BRCA = read.csv('Data4Analysis/Breast_clinicalData_07.07.21.tsv', sep = '\t')
Paths = read.csv('~/Documents/MSKCC/05_IMPACT40K/Data/Signed_out/Facets_annotated.cohort.txt', sep = '\t')

random.Normals = BRCA[sample(nrow(BRCA), size = 1020, replace = F), 'Sample.ID']
paths.Normals = Paths[which(Paths$tumor_sample %in% random.Normals), 'counts_file']
write.table(paths.Normals, file = 'DataProcessed/PON_BRCA_Paths.txt', col.names = F, row.names = F, quote = F)


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