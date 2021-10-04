## Data selection; BRCA cohort; benchmark DryClean;
## 
## 09/18/2021
## chris-kreitzer

rm(list = ls())
set.seed(136)


## Input 
cohort_compendium = read.csv('~/Documents/MSKCC/04_Y_chromo_loss/Data/WES_cBio_ID.match.txt', sep = '\t')
BRCA = cohort_compendium[which(cohort_compendium$Tumor_Type == 'Breast Invasive Ductal Carcinoma'), ]
PP_annotations = read.csv('~/Documents/MSKCC/05_IMPACT40K/Data/Signed_out/Facets_annotated.cohort.txt', sep = '\t')

## Processing
BRCA_annotation = PP_annotations[which(PP_annotations$tumor_sample %in% BRCA$DMP_Sample_ID), ]
BRCA_annotation = merge(BRCA_annotation, BRCA, by.x = 'tumor_sample', by.y = 'DMP_Sample_ID', all.x = T)

#' high purity samples
high_purity.cohort = BRCA_annotation[which(BRCA_annotation$purity > 0.65), ]
low_purity.cohort = BRCA_annotation[which(BRCA_annotation$purity < 0.35), ]


preRun_BRCA.cohort_full = rbind(high_purity.cohort, low_purity.cohort)
preRun_BRCA.cohort_short = preRun_BRCA.cohort_full[, c('tumor_sample', 'counts_file', 'path', 'purity', 'wgd', 
                                                       'ploidy', 'n_amps', 'n_homdels', 'frac_homdels',
                                                       'CMO_Sample_ID', 'Tumor_Type', 'Facet_Path', 'Facet_Countfile')]

BRCA_workingCohort_MSK = list(full_data = preRun_BRCA.cohort_full,
                              short_data = preRun_BRCA.cohort_short)

saveRDS(BRCA_workingCohort_MSK, file = '~/Documents/GitHub/DryClean_Facets/Data_out/BRCA_workingCohort_MSK.rds')



# write.table(x = data.frame(sample = preRun_BRCA.cohort_short$counts_file), 
#             file = '~/Documents/GitHub/DryClean_Facets/Data_out/Tumor_Paths.txt', sep = '\t', row.names = F)

#' out


#' 09/19/2021
#' update the cohort - based on conversion
old_cohort = readRDS('Data_out/BRCA_workingCohort_MSK.rds')
IDs = old_cohort$short_data$tumor_sample
successfully_converted = list.files(path = 'Tumor_samples/', pattern = '^counts.*', full.names = F, all.files = T)
successfully_converted = substr(successfully_converted, start = 17, stop = 33)

successfully_converted = old_cohort$short_data[which(old_cohort$short_data$tumor_sample %in% successfully_converted), ]
new_cohort = append(old_cohort, list(successfully_converted))
names(new_cohort)[3] = 'data'
saveRDS(successfully_converted, file = 'Data_out/BRCA_workingCohort_MSK.rds')

#' fetch normal samples
cohort = readRDS('Data_out/BRCA_workingCohort_MSK.rds')
normal_counts = data.frame(sample = cohort$counts_file)
write.table(normal_counts, file = 'Normal_samples/Normals.txt', sep = '\t', quote = F, row.names = F)
