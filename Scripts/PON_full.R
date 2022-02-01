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

