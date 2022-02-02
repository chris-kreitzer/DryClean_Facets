## Create a Panel of Normals (PON)
##
## Learn, whether different tools provide us with different count matrices
## Learn, about the noise in the input data
## Learn, whether a reduced AND/OR binned PON make a difference
##
## 02/01/2022
## chris-kreitzer


clean()
setup(working.path = '~/Documents/GitHub/DryClean_Facets/')


## Libraries and Dependencies:
devtools::install_github('mskilab/dryclean')
# devtools::install_github('mskilab/gUtils')
devtools::install_github("mskilab/bamUtils")
BiocManager::install('S4Vectors')
BiocManager::install('GenomicAlignments')
library(S4Vectors)
library(gUtils)
library(dryclean)
library(tidyverse)
library(pbmcapply)
library(data.table)


## Check the average Coverage across a panel of normal samples
## here I just used a random set of 60 IMPACT samples
NCOV = c()
for(i in list.files('~/Desktop/mnt/ATMcountdata/', full.names = T)){
  input = vroom::vroom(i)
  ii = input$File1R + input$File1A
  ii = ii[which(ii > 30)]
  ii_mean = mean(ii)
  NCOV = c(NCOV, ii_mean)
}




