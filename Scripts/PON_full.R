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
library(GenomicRanges)

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
erbb2 = input[which(input$Chromosome == 17 & 
                      input$Position >= gene.start & 
                      input$Position <= gene.end), 
              c('Chromosome', 'Position', 'NOR.DP', 'TUM.DP')]

## ERBB2 coverage from normal sample: Marcin Imiliensky approach
#' general annotations
genomic_features = import("gencode.v19.chr_patch_hapl_scaff.annotation.gtf.gz")
exons = genomic_features %Q% (type == 'exon')
genes = genomic_features %Q% (type == 'gene')
gt.ge = track.gencode(genes = 'ERBB2', 
                      bg.col = 'white',
                      cds.col = alpha('red', 0.9),
                      utr.col = 'black', 
                      cex.label = 1,
                      st.col = 'orange',
                      en.col = 'grey55')

erbb2_GR = makeGRangesFromDataFrame(df = erbb2,
                                    keep.extra.columns = T,
                                    start.field = 'Position',
                                    end.field = 'Position',
                                    ignore.strand = T)
names(erbb2_GR) = NULL
dc.dcb = gTrack(data = erbb2_GR, 
                y.field = 'NOR.DP', 
                ygap = 0.9, 
                col = 'black', 
                name = 'Normal Coverage',
                circles = F, 
                lwd.border = 2, 
                y1 = 1200,
                xaxis.width = 1, 
                height = 20, yaxis.pretty = T, 
                yaxis.cex = 1, formatting = T)

win = (genes %Q% (gene_name == 'ERBB2') + 1e2) %&% exons %Q% (1)
plot(c(gt.ge, dc.dcb), win, col = 'black', border = 1)



###############################################################################
## Input data: fetched from juno
BRCA_PON_list = readRDS('~/Documents/MSKCC/07_FacetsReview/DataProcessed/BRCA_PON_list.rds')
BRCA_PON_df = data.table::rbindlist(BRCA_PON_list)

#' automate marker selection for proper dimensions in PON
#' Note, that n (marker-bins) x m(samples) need to be equal among all normal samples
max_PON = function(normal_samples){
  
  #' with this function we create a max representation of our Normals
  
  max_n = max(table(normal_samples$sample))
  min_n = min(table(normal_samples$sample))
  
  input_list = normal_samples
  input_list$duplication = paste(input_list$Chromosome, input_list$Position, sep = ';')
  bins_PON = data.frame()
  sample_PON = data.frame()
  message('Creating a maximal representation of normals (PON)')
  
  #' search for samples which have duplicated entries (chromosome*position) and discard them
  #' write a little function to catch samples with duplicated entries
  dupli_events = function(data, sample){
    if(any(duplicated(data$duplication[which(data$sample == sample)]))){
      unique(data$sample[which(data$sample == sample)])
    }
  }
  
  x = sapply(unique(input_list$sample), function(x) dupli_events(data = input_list, sample = x))
  x_reduced = Filter(Negate(is.null), x)
  x_reduced = as.character(unlist(x_reduced))
  
  #' subset input list; to remove samples with duplicated entries
  input_list = input_list[!input_list$sample %in% x_reduced, ]
  rm(x, x_reduced)
  message('Sample Quality Control ended')
  
  
  #' loop through list and select common positions
  n.PON = length(unique(input_list$sample))
  max_bins = input_list$sample[which.max(table(input_list$sample))]
  max_bins_sample = input_list[input_list %in% max_bins, ]
  matrix.table.keep = data.frame(loc = max_bins_sample$duplication)
  
  #' sample-wise listing of postions
  locations.out = data.frame()
  for(patient in unique(input_list$sample)){
    print(patient)
    if(all(matrix.table.keep$loc %in% input_list$duplication[which(input_list$sample == patient)])){
      table.out = input_list[which(input_list$sample == patient & input_list$duplication %in% matrix.table.keep$loc), ]
    } else {
      table.out = input_list[which(input_list$sample == patient & input_list$duplication %in% matrix.table.keep$loc), ]
      missing = setdiff(matrix.table.keep$loc, input_list$duplication[which(input_list$sample == patient)])
      missing.df = data.frame(duplication = missing,
                              sample = patient)
      missing.df = separate(missing.df, 
                            col = duplication,
                            into = c('Chromosome', 'Position'),
                            sep = ';',
                            remove = F)
      
      #' add artificial data for missing positions; in this case just 1
      missing.df$NOR.DP = 1
      missing.df$NOR.RD = 1
      table.out = rbind(table.out, missing.df)
    }
    locations.out = rbind(locations.out, table.out)
  }
  
  
  #' prepare the final output
  PON_out = as.data.frame(do.call('cbind', split(locations.out[, c('NOR.DP')], locations.out$sample)))
  row.names(PON_out) = matrix.table.keep$loc
  
  #' mean-normalization
  mean_normalization = function(x){
    x / mean(x)
  }
  
  message('Starting mean-normalization')
  
  PON_normalized = apply(PON_out, 2, mean_normalization)
  
  #' return object
  return(list(selcted_bins = matrix.table.keep$loc,
              PON_normalized = PON_normalized))
}





