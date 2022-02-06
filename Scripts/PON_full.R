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
library(rtracklayer)


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

erbb2 = input[which(input$Chromosome == 17 & input$Position >= gene.start & input$Position <= gene.end), 
              c('Position', 'NOR.DP', 'TUM.DP') ]



###############################################################################
## Input data: fetched from juno
BRCA_PON_list = readRDS('~/Documents/MSKCC/07_FacetsReview/DataProcessed/BRCA_PON_list.rds')
BRCA_PON_df = data.table::rbindlist(BRCA_PON_list)

#' automate marker selection for proper dimensions in PON
#' Note, that n (marker-bins) x m(samples) need to be equal among all normal samples



library(rtracklayer)
this.gr = import("gencode.v19.chr_patch_hapl_scaff.annotation.gtf.gz")
metadata = data.frame(sampleid = c('P-0000584-T03-IM6', 'P-0003195-T02-IM6'),
                      cov = c('/Users/chriskreitzer/Documents/GitHub/DryClean_Facets/Normal_samples/countsMerged____P-0000584-T03-IM6_P-0000584-N01-IM6.dat.gz.rds',
                              '/Users/chriskreitzer/Documents/GitHub/DryClean_Facets/Normal_samples/countsMerged____P-0003195-T02-IM6_P-0003195-N01-IM6.dat.gz.rds'),
                      dryclean = c('/Users/chriskreitzer/Documents/GitHub/DryClean_Facets/Tumor_cleaned/P-0000584-T03-IM6_P-0000584-N01-IM6_drycleaned.rds',
                                   '/Users/chriskreitzer/Documents/GitHub/DryClean_Facets/Tumor_cleaned/P-0003195-T02-IM6_P-0003195-N01-IM6_drycleaned.rds'))


cov = metadata[1, 'cov'] %>% readRDS()
dc = metadata[1, 'dryclean'] %>% readRDS()
gt.dcb = gTrack(dc, 'background', circle=TRUE, lwd.border=0.8)

exons = this.gr %Q% (type == 'exon')
genes = this.gr %Q% (type == 'gene')
gt.ge = track.gencode()
gtr = gTrack(reduce(cov))
win = (genes %Q% (gene_name == 'ERBB2') + 1e2) %&% exons %Q% (1)
plot(c(gt.ge, gt.dcb, gtr), win, col = 'magenta' , border = '60')

gr <- GRanges(seqnames = Rle(c("chr1" , "chr2" , "chr1" , "chr3") ,
                             c(1,3,2,4)), ranges = IRanges(c(1,3,5,7,9,11,13,15,17,19) ,
                                                           end = c(2,4,6,8,10,12,14,16,18,20),
                                                           names = head(letters,10)),
              GC=seq(1,10,length=10),
              name=seq(5,10,length=10))

graph = matrix(0 , nrow = 10 , ncol = 10)
graph[1,3]=1
graph[1,10]=1
graph[2,5]=1
graph[2,8]=1
graph[3,5]=1
graph[4,1]=1
graph[4,2]=1
graph[4,6]=1
graph[4,9]=1
graph[5,1]=1
graph[5,2]=1
graph[5,4]=1
graph[8,1]=1
graph[8,2]=1
graph[9,1]=1
graph[10,1]=1

plot(gTrack(gr , edges = graph , stack.gap = 5),
     colormaps = NULL, # (named) list same length as data
     height = 50,
     ygap = 0.2,
     bg.col = 'white',
     ylab = 'hallo',
     stack.gap = 0,
     cex.label = 2,
     gr.cex.label = 2 *0.8,
     gr.srt.label = 0,
     col = NA,
     border = NA,
     angle = 15,
     name = "hallo",
     gr.colorfield = NA,
     y.quantile = 0.01, ## if y0 or y1 is not specified then will draw between y.quantile and 1-y.quantile of data
     y.cap = T, ## whether to cap values at y0, y1 (only relevant if y.field specified)
     lwd.border = 1,
     hadj.label = 1,
     vadj.label = 0.5,
     smooth = NA, ## smooth with running mean with this window
     round = NA, ## round the output of running mean to this precision
     ywid = NA,
     ypad = 0,
     seqinfo = NA,
     circles = TRUE,
     lines = TRUE,
     bars = FALSE,
     draw.paths = FALSE,
     path.col = 'green',
     path.lwd = 3 )





