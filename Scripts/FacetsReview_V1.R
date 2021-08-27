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
set.seed(111)

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
median.seg1 = median(chr1$cnlr[which(chr1$maploc < 105977093)])
median.seg2 = median(chr1$cnlr[which(chr1$maploc >= 105977093)])

chr1$plot.dummy = seq(1, nrow(chr1), 1)
chr1$segclust = as.factor(as.character(as.numeric(chr1$segclust)))
logR = ggplot(chr1, aes(x = plot.dummy, y = cnlr, color = segclust, alpha = as.character(het))) + 
        geom_point() +
        scale_color_manual(values = c('5' = 'blue', '9' = 'red')) +
        scale_y_continuous(expand = c(0, 0),
                           limits = c(-2, 2),
                           breaks = seq(-2, 2, 0.5)) +
        scale_x_discrete(expand = c(0.01,0.01)) +
        annotate('segment', x = 0, xend = chr1$plot.dummy[which.max(chr1$seg0)] - 1,
                 y = median.seg1, yend = median.seg1, colour = 'blue') +
        annotate('segment', x = chr1$plot.dummy[which.max(chr1$seg0)] - 1, 
                 xend = max(chr1$plot.dummy),
                 y = median.seg2, yend = median.seg2, colour = 'red') +
        theme_bw() +
        theme(axis.text.x = element_blank(),
              axis.text.y = element_text(size = 12, face = 'bold', colour = 'black'),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()) +
        geom_hline(yintercept = 0, color = 'grey', linetype = 'dashed') +
        labs(x = '', y = 'Log Ratio / LogR')



logR

#' plot the logOR of chromosome 1
baf = facets_fit$cncf[which(facets_fit$cncf$chrom == 1), ]
mafR_1 = baf$mafR[which(baf$chrom == 1 & which.min(baf$segclust))]
mafR_a = mafR_1[1]
mafR_b = mafR_1[2]

logOR = ggplot(chr1, aes(x = plot.dummy, y = valor, color = segclust)) + 
        geom_point() +
        scale_color_manual(values = c('5' = 'blue', '9' = 'red')) +
        annotate('segment', x = 0, xend = chr1$plot.dummy[which.max(chr1$seg0)] - 1,
                 y = mafR_a, yend = mafR_a, colour = 'blue') +
        annotate('segment', x = 0, xend = chr1$plot.dummy[which.max(chr1$seg0)] - 1,
                 y = -mafR_a, yend = -mafR_a, colour = 'blue') +
        annotate('segment', x = chr1$plot.dummy[which.max(chr1$seg0)] - 1, xend = which.max(chr1$plot.dummy),
                 y = mafR_b, yend = mafR_b, colour = 'red') +
        annotate('segment', x = chr1$plot.dummy[which.max(chr1$seg0)] - 1, xend = which.max(chr1$plot.dummy),
                 y = -mafR_b, yend = -mafR_b, colour = 'red') +
        scale_y_continuous(expand = c(0, 0),
                           limits = c(-2, 2),
                           breaks = seq(-2, 2, 0.5)) +
        scale_x_discrete(expand = c(0.01,0.01)) +
        
        theme_bw() +
        theme(axis.text.x = element_blank(),
              axis.text.y = element_text(size = 12, face = 'bold', colour = 'black'),
                                         panel.grid.major = element_blank(),
                                         panel.grid.minor = element_blank()) +
        geom_hline(yintercept = 0, color = 'grey', linetype = 'dashed') +
        labs(x = 'genomic coordinates', y = 'log-ODDS-ratio / logOR')

logR / logOR  
        

#' make a modified version of the plot; for easier handling and interpretation
#' how does a diploid genome look like?
logR_modi = ggplot(baf, aes(x = seg, y = cnlr.median, color = as.character(seg), 
                fill = as.character(seg))) +
        geom_bar(stat = 'identity') +
        scale_color_manual(values = c('1' = 'blue', '2' = 'red')) +
        scale_fill_manual(values = c('1' = 'blue', '2' = 'red')) +
        geom_hline(yintercept = 0, color = 'grey', linetype = 'dashed') +
        scale_y_continuous(expand = c(0, 0),
                           limits = c(-2, 2),
                           breaks = seq(-2, 2, 0.5)) +
        scale_x_discrete(expand = c(0.01,0.01)) +
        
        theme_bw() +
        theme(axis.text.x = element_blank(),
              axis.text.y = element_text(size = 12, face = 'bold', colour = 'black'),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()) +
        labs(x = '', y = 'Log-Ratio / abs(deviation) from zero')


#' make the corresponding BAF plot     
logOR_modi = ggplot(baf, aes(x = seg, y = abs(mafR), color = as.character(seg), 
                fill = as.character(seg))) +
        geom_bar(stat = 'identity') +
        scale_color_manual(values = c('1' = 'blue', '2' = 'red')) +
        scale_fill_manual(values = c('1' = 'blue', '2' = 'red')) +
        geom_hline(yintercept = 0, color = 'grey', linetype = 'dashed') +
        scale_y_continuous(expand = c(0, 0),
                           limits = c(0, 2),
                           breaks = seq(0, 2, 0.5)) +
        scale_x_discrete(expand = c(0.01,0.01)) +
        
        theme_bw() +
        theme(axis.text.x = element_blank(),
              axis.text.y = element_text(size = 12, face = 'bold', colour = 'black'),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()) +
        
        # annotate('rect', xmin = 0, xmax = nrow(baf), ymin = 0, ymax = 0.2,
        #          color = 'grey88', fill = 'grey88', alpha = 0.2) +
        labs(x = '', y = 'abs(deviation) logOR ratio')

logR_modi / logOR_modi







