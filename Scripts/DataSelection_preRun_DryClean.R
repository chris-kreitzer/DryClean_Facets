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
high_purity.cohort = BRCA_annotation[which(BRCA_annotation$purity > 0.70), ]
low_purity.cohort = BRCA_annotation[which(BRCA_annotation$purity < 0.30), ]

