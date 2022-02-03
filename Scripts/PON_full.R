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
.rs.restartR()
org.Hs.egACCNUM

## Plot exon structure of ERBB2
org.Hs.eg.db

genome <- TxDb.Hsapiens.UCSC.hg19.knownGene
genic.regions <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
# the plec gene
erbb2_gene = genes(genome)[which(genes(genome)$gene_id == 2064),]

# get the exons with the gene coordinates
plec_exons = subsetByOverlaps(exons(genome), erbb2_gene)

# same for transcripts
plec_t = subsetByOverlaps(transcripts(genome), plec_gene)



#store the first six keys
my_keys = head(keys(org.Hs.eg.db))

keytypes(org.Hs.eg.db)

columns(org.Hs.eg.db)

#selecting
select(org.Hs.eg.db,
       keys = my_keys,
       columns=c("ENTREZID","SYMBOL","GENENAME"),
       keytype="ENTREZID")




ERBB2_sub = input[which(input$Chromosome == 17 & input$Position > 37854492 & input$Position < 37886297), ]




ERBB2_sub
plot(ERBB2_sub$TUM.DP)

exon1 = ERBB2_sub[which(ERBB2_sub$Position > 37856292 & ERBB2_sub$Position < 37856764), ]
plot(exon1$TUM.DP)
plot(exon1$TUM.DP ~ exon1$Position)



library(ggplot2)
ERBB2_coord$st = seq(1, nrow(ERBB2_coord)*2, 2)
ERBB2_coord$stop = seq(2, nrow(ERBB2_coord)*2, 2)



ggplot(ERBB2_coord, aes(x = st, y = st, width = start - end)) +
  geom_tile()









