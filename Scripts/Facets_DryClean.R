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


##-----------------------------------------------------------------------------
## Look into the Facets calls;
Facets_Dryclean = list.files(path = '~/Documents/MSKCC/07_FacetsReview/DryClean/Facets_DryClean/', full.names = T)
Facets_Clean = data.frame()
for(i in 1:length(Facets_Dryclean)){
  print(i)
  data_in = readRDS(Facets_Dryclean[i])
  purity = data_in$Facets$purity
  ploidy = data_in$Facets$ploidy
  diplogr = data_in$Facets$dipLogR
  flags = data_in$Facets$em_flags
  if(is.null(flags)) flags = NA
  if(length(flags) > 1) flags = paste(flags, collapse = ';')
  segments = nrow(data_in$Facets$segs)
  wgd = data_in$QC$wgd
  fga = data_in$QC$fga
  qc = data_in$QC$facets_qc
  
  out = data.frame(sample = basename(Facets_Dryclean[i]),
                   purity = purity,
                   ploidy = ploidy,
                   diplogr = diplogr,
                   flags = flags,
                   segments = segments,
                   wgd = wgd,
                   fga = fga,
                   qc = qc)
  Facets_Clean = rbind(Facets_Clean, out)
  
}


##-----------------------------------------------------------------------------
## Facets Original
Facets_original = list.files(path = '~/Documents/MSKCC/07_FacetsReview/DryClean/Tumor_Facets/', full.names = T)
Facets_Original = data.frame()
for(i in 1:length(Facets_original)){
  print(i)
  data_in = readRDS(Facets_original[i])
  purity = data_in$facets_output$purity
  ploidy = data_in$facets_output$ploidy
  diplogr = data_in$facets_output$dipLogR
  flags = data_in$facets_output$em_flags
  if(is.null(flags)) flags = NA
  if(length(flags) > 1) flags = paste(flags, collapse = ';')
  segments = nrow(data_in$facets_output$segs)
  wgd = data_in$qc$wgd
  fga = data_in$qc$fga
  qc = data_in$qc$facets_qc
  
  out = data.frame(sample = basename(Facets_original[i]),
                   purity = purity,
                   ploidy = ploidy,
                   diplogr = diplogr,
                   flags = flags,
                   segments = segments,
                   wgd = wgd,
                   fga = fga,
                   qc = qc)
  
  Facets_Original = rbind(Facets_Original, out)
  
}

  
a = readRDS(Facets_original[1])
a$facets_output  