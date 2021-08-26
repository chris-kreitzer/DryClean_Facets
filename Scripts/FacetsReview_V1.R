## Demonstrative Purposes of Facets; allele-specific alteration pattern:
## Take a sample where we have IMPACT Targeted-Panel and WES data and a relatively pure sample
## 
## Start: 06/21/2021
## chris kreitzer


#' sample selection:
#' see description on word-file
#' P-0020818-T01-IM6

rm(list = ls())
setwd('~/Documents/GitHub/DryClean_Facets/')


## Libraries:
library(pctGCdata)
library(patchwork)


## Input:
FacetsAnnotation = read.csv('~/Documents/MSKCC/05_IMPACT40K/Data/Signed_out/Facets_annotated.cohort.txt', sep = '\t')
TP_sample = read.csv('Data4Analysis/countsMerged____P-0020818-T01-IM6_P-0020818-N01-IM6.dat.gz')

#' /juno/work/ccs/shared/resources/impact/facets/all/P-00208/P-0020818-T01-IM6_P-0020818-N01-IM6/	
#' /juno/work/tempo/wes_repo/Results/v1.3.x/cohort_level/MSKWESRP/somatic/s_C_001451_P001_d__s_C_001451_N001_d/facets/s_C_001451_P001_d__s_C_001451_N001_d/s_C_001451_P001_d__s_C_001451_N001_d.snp_pileup.gz

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Running Facets default parameters on described GBM sample yields the following picture:

TP_sample_facets = facets::readSnpMatrix(filename = '~/Documents/MSKCC/07_FacetsReview/Data4Analysis/countsMerged____P-0020818-T01-IM6_P-0020818-N01-IM6.dat.gz')
pre_process = facets::preProcSample(rcmat = TP_sample_facets, 
                                    ndepth = 35, 
                                    het.thresh = 0.25, # vaf threshold for het-SNP in normal sample)
                                    snp.nbhd = 250, # insert size; 250 suggested for IMPACT TP  
                                    cval = 25,
                                    gbuild = 'hg19',
                                    hetscale = T)
post_process = facets::procSample(x = pre_process,
                                  cval = 150, 
                                  min.nhet = 15, # min number of het-SNPs per segment required for segmentation; Hotelling's)
                                  dipLogR = NULL) 

facets_fit = facets::emcncf(x = post_process, 
                            trace = TRUE, 
                            maxiter = 1e6)




##-- Question 1: Heterozygous loss at chromosome 1
#' extract information from Chromosome 1 and look at the logOR from Facets and BAF and the relation
chr1 = post_process$jointseg[which(post_process$jointseg$chrom == 1), ]
median.seg1 = median(chr1$cnlr)


View(chr1)        
dim(chr1)
max(chr1$cnlr, na.rm = T)
min(chr1$cnlr, na.rm = T)


a = 
min(a$cnlr)
max(a$cnlr)

View(a)

chr1_raw_seg1 = chr1_raw[which(chr1_raw$seg == 1), ]
chr1_raw_hetero = chr1_raw_seg1[!is.na(chr1_raw_seg1$valor), ]

#' get the heterozygous positions:
chr1_heterozygous = TP_sample[which(TP_sample$Chromosome == 1 & TP_sample$Position %in% chr1_raw_hetero$maploc), ]
chr1_heterozygous$BAF = chr1_heterozygous$File2A / (chr1_heterozygous$File2A + chr1_heterozygous$File2R)




pre_process$pmat
pre_process$seg.tree
pre_process$jointseg









