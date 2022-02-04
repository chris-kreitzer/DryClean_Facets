## select 1,000 random Normal Samples where the PON
## will be created.
## this script will be run on the terminal (juno)
## BRCA data (gene panel 468) was fetched from cBIO: 07/07/2021


BRCA = read.csv('Data4Analysis/Breast_clinicalData_07.07.21.tsv', sep = '\t')
Paths = read.csv('~/Documents/MSKCC/05_IMPACT40K/Data/Signed_out/Facets_annotated.cohort.txt', sep = '\t')

random.Normals = BRCA[sample(nrow(BRCA), size = 1020, replace = F), 'Sample.ID']
paths.Normals = Paths[which(Paths$tumor_sample %in% random.Normals), 'counts_file']
write.table(paths.Normals, file = 'DataProcessed/PON_BRCA_Paths.txt', col.names = F, row.names = F, quote = F)

#' fetch coordinates from Normal samples; 
#' this script will be submitted to LFS on the juno-cluster
library(facets)
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