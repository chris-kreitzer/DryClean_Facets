library(FacetsDC)


#' create table for FacetsDC
paths = read.csv('~/Documents/MSKCC/07_FacetsReview/DryClean/DryClean_Facets_table.txt', sep = '\t')
paths$count_file = paste0('~/Documents/MSKCC/07_FacetsReview/DryClean/Tumor_Countfiles/', basename(paths$count_file))
decomposed = list.files(path = '~/Documents/MSKCC/07_FacetsReview/DryClean/Tumor_Decomposed/', full.names = T)
decomposed = data.frame(value = decomposed, 
                        original = substr(x = basename(decomposed), 
                                          start = 1, 
                                          stop = nchar(basename(decomposed)) - 8))

paths = merge(paths, decomposed, by.x = 'sample', by.y = 'original', all.x = T)
paths$tumor_clean = NULL
colnames(paths)[5] = 'tumor_clean'
write.table(paths, file = '~/Documents/MSKCC/07_FacetsReview/DryClean/DryClean_Facets_table.txt', sep = '\t', quote = F, row.names = F)

##-----------------------------------------------------------------------------
## run FacetsDC on decomposed samples and count files
clean()
gc()
source('Scripts/FacetsQC.R')
library(FacetsDC)
samples = read.csv('~/Documents/MSKCC/07_FacetsReview/DryClean/DryClean_Facets_table.txt', sep = '\t')

for(i in 1:nrow(samples)){
  print(i)
  out = FacetsDC::run_facets_cleaned(read_counts = samples$count_file[i],
                                     read_cleaned = samples$tumor_clean[i],
                                     MODE = 'union', 
                                     cval = 150)
  QC = facets_fit_qc(facets_output = out)
  return_object = list(Facets = out, QC = QC)
  saveRDS(return_object, file = paste0('~/Documents/MSKCC/07_FacetsReview/DryClean/Facets_DryClean/', samples$sample[i], '.rds'))
}


