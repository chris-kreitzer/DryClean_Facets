#' ## PACKAGE Testing:
#' 
#' rm(list = ls())
#' .rs.restartR()
set.seed(100)
setwd('~/Documents/GitHub/DryClean_Facets/Scripts/DryClean_Segmentation/facets-0.6.1/')

#' testing the updated package
install.packages('~/Documents/GitHub/DryClean_Facets/Scripts/DryClean_Segmentation/facets-0.6.1/',
                 type = "source", repos = NULL)

library(FacetsDC)

