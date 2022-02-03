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
# devtools::install_github("mskilab/bamUtils")
# BiocManager::install('S4Vectors')
# BiocManager::install('GenomicAlignments')
# remotes::install_github("mskcc/facets", build_vignettes = TRUE)
devtools::install_github("mskilab/gTrack")
devtools::install_github("mskilab/skidb")
devtools::install_github("mskilab/skitools")


library(S4Vectors)
library(gUtils)
library(dryclean)
library(tidyverse)
library(pbmcapply)
library(data.table)
library(facets)
library(pctGCdata)
library("GenomicFeatures")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
# BiocManager::install("org.Hs.eg.db")
library("org.Hs.eg.db")
# BiocManager::install('Repitools')
library(Repitools)
library(patchwork)
library('skitools')
library('gTrack')
library('skidb')


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

## Check the sequencing distribution of ERBB2 
ERBB2_coord = read.csv(file = 'Data_out/ERBB2_Probes.txt', sep = '\t')
input = facets::readSnpMatrix('~/Documents/MSKCC/07_FacetsReview/Tumor_countsFile/countsMerged____P-0000584-T03-IM6_P-0000584-N01-IM6.dat.gz')

## Plot exon structure of ERBB2
x = org.Hs.egACCNUM
mapped_genes = mappedkeys(x)
ID_match = select(org.Hs.eg.db,
       keys = mapped_genes,
       columns = c("ENTREZID","SYMBOL","GENENAME"),
       keytype = "ENTREZID")
xx = as.list(org.Hs.egALIAS2EG)


genome = TxDb.Hsapiens.UCSC.hg19.knownGene
genic.regions = genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
erbb2_gene = genes(genome)[which(genes(genome)$gene_id == 2064), ]

# get the exons with the gene coordinates
erbb2_exons = subsetByOverlaps(exons(genome), erbb2_gene)
erbb2_exons = annoGR2DF(erbb2_exons) 
