## DryClean; Creating Panel of Normal:
## 
## Here, we will start focusing on Breast adenocarcinoma samples
## We will be using IMPACT gene panel 468 for this appraoch;
## Our PON will consist of 100 normal samples
## 
## Full detail: Data queried at 07/07/2021 from cBIO
## Breast Cancer; Breast Invasive Ductal Carcinoma; Gene Panel 468 (n = 2,705 samples)
## 
## 07/07/2021
## chris kreitzer


rm(list = ls())
.rs.restartR()
setwd('~/Documents/MSKCC/07_FacetsReview/')
set.seed(111)


## Libraries and Dependencies
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = TRUE)
# devtools::install_github('mskilab/dryclean')
# devtools::install_github('mskilab/gUtils')
# devtools::install_github('mskilab/fragCounter')
# devtools::install_github("mskilab/bamUtils")
# devtools::install_github('mskilab/fragCounter')   
# BiocManager::install('S4Vectors')
# BiocManager::install('GenomicAlignments')
# library(S4Vectors)
library(gUtils)
library(dryclean)


## Input data:
BRCA = read.csv('Data4Analysis/Breast_clinicalData_07.07.21.tsv', sep = '\t')
Paths = read.csv('~/Documents/MSKCC/05_IMPACT40K/Data/Signed_out/Facets_annotated.cohort.txt', sep = '\t')
BRCA_PON_list = readRDS('DataProcessed/BRCA_PON_list.rds')


## Data wrangling and processing
#' create a Panel of Normal
random.Normals = BRCA[sample(nrow(BRCA), size = 100, replace = F), 'Sample.ID']
paths.Normals = Paths[which(Paths$tumor_sample %in% random.Normals), 'counts_file']
write.table(paths.Normals, file = 'DataProcessed/PON_BRCA_Paths.txt', col.names = F, row.names = F, quote = F)

#' fetch coordinates from Normal samples; 
#' this script will be submitted to LFS on the juno-cluster
# library(facets)
# PON = read.csv('/juno/home/kreitzec/DryClean/PON_BRCA_Paths.txt', sep = '\t', header = F, row.names = F)
# 
# BRCA_PON_list = list()
# for(i in unique(PON)){
#   data.in = facets::readSnpMatrix(i)
#   data.processed = data.in[which(data.in$NOR.DP >= 35), c('Chromosome', 'Position', 'NOR.DP', 'NOR.RD')]
#   data.processed$sample = substr(x = basename(i), start = 17, stop = 33)
#   BRCA_PON_list[[i]] = data.processed
# }
# 
# saveRDS(object = BRCA_PON_list, file = '/juno/home/kreitze/DryClean/BRCA_PON_list.rds')


#' fetch all the positions which are shared among all n = 99 normal samples


#' make the frequency table per chromosome and consider also
#' Positions with < 99 (90) - meaning that 90% must contain postion




a = lapply(BRCA_PON_list, function(x) Reduce(intersect, list(x$Position)))
library(magrittr)
i = lapply(BRCA_PON_list, function(i) i[,c(1, 2)]) %>% do.call(as.data.frame(xtabs(~1 + 2)))



result = lapply(BRCA_PON_list, "[", , c('Position', 'Chromosome'))

result[[1]]
head(vec)
result[[1]]
head(str(result))
result[[1]]

all = c()
for(i in seq_along(result)){
  u = unlist(result[[i]])
  all = c(all, u)
}



BRCA_PON_list[[1]]
u = result[[1]]
head(u)
y = as.data.frame(xtabs(~u$Chromosome + u$Position))
head(y)
dim(y)
View(y)















