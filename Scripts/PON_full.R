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
# devtools::install_github('mskilab/dryclean')
# devtools::install_github('mskilab/gUtils')
# devtools::install_github("mskilab/bamUtils")
# BiocManager::install('S4Vectors')
# BiocManager::install('GenomicAlignments')
# remotes::install_github("mskcc/facets", build_vignettes = TRUE)
# devtools::install_github("mskilab/gTrack")
# devtools::install_github("mskilab/skidb")
# devtools::install_github("mskilab/skitools")
# BiocManager::install("org.Hs.eg.db")
# BiocManager::install('Repitools')

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
library(dryclean)
library("org.Hs.eg.db")
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
ERBB2_probes = read.csv(file = 'Data_out/ERBB2_Probes.txt', sep = '\t')
input = facets::readSnpMatrix('~/Documents/MSKCC/07_FacetsReview/Tumor_countsFile/countsMerged____P-0000584-T03-IM6_P-0000584-N01-IM6.dat.gz')


## Extract genomic features from a TxDb-like object
## Plot exon structure of ERBB2
ii = org.Hs.egACCNUM
mapped_genes = mappedkeys(ii)
ID_match = select(org.Hs.eg.db,
       keys = mapped_genes,
       columns = c("ENTREZID","SYMBOL","GENENAME"),
       keytype = "ENTREZID")

#' retrieve exons from specific gene ERBB2
#' Columns to select: 
columns(TxDb.Hsapiens.UCSC.hg19.knownGene)
erbb2_exons = exons(TxDb.Hsapiens.UCSC.hg19.knownGene,
                    columns = c("EXONSTART", 'EXONEND', 'GENEID'),
                    filter = list(gene_id = ID_match$ENTREZID[which(ID_match$SYMBOL == 'ERBB2')]))

erbb2_exons = annoGR2DF(erbb2_exons)

#' gene start
gene.start = 37844347
gene.end = 37884911

erbb2 = input[which(input$chrom == 17, input$Postion)]


head(erbb2_exons)
head(ERBB2_probes)

ggplot() +
  geom_rect(aes(xmin = 37856492, xmax = 37856564, ymin = -2, ymax = 0), col = 'blue', fill = 'blue') +
  geom_rect(aes(xmin = 37844393, xmax = 37844531, ymin = 0.5, ymax = 1))









###############################################################################
## Input data:
BRCA = read.csv('Data4Analysis/Breast_clinicalData_07.07.21.tsv', sep = '\t')
Paths = read.csv('~/Documents/MSKCC/05_IMPACT40K/Data/Signed_out/Facets_annotated.cohort.txt', sep = '\t')
BRCA_PON_list = readRDS('DataProcessed/BRCA_PON_list.rds')


## Data wrangling and processing
#' create a Panel of Normal; n = 1,000
random.Normals = BRCA[sample(nrow(BRCA), size = 1020, replace = F), 'Sample.ID']
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


BRCA_PON_list = readRDS('DataProcessed/BRCA_PON_list.rds')

#' replace with Rbindlist
BRCA_PON_df = rbindlist(BRCA_PON_list)

#' automate marker selection for proper dimensions in PON
#' Note, that n (marker-bins) x m(samples) need to be equal among all normal samples
library(EnsDb.Hsapiens.v86)
exons(x, ...)
## S4 method for signature 'TxDb'
exons(x, columns="exon_id", filter=NULL, use.names=FALSE)
BiocManager::install("EnsDb.Hsapiens.v86")


gene_id

exons(TxDb.Hsapiens.UCSC.hg19.knownGene, columns = c("EXONID", "TXNAME"),
      filter=list(gene_id=2064))



a = data.frame(x = seq(1, 50, 1),
               y = rnorm(n = 50, mean = 10, sd = 2))

ggplot(a) +
  geom_rect(aes(xmin = 1, xmax = 50, ymin = -2, ymax = 0), col = 'blue', fill = 'blue') +
  geom_rect(aes(xmin = 20, xmax = 30, ymin = 0.5, ymax = 1)) +
  geom_point(aes(x = x, y = y))


head(a)

            