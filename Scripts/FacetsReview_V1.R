## Demonstrative Purposes of Facets; allele-specific alteration pattern:
## Take a sample where we have IMPACT Targeted-Panel and WES data and a relatively pure sample
## 
## Start: 06/21/2021
## chris kreitzer


#' sample selection:
#' see description on word-file
#' P-0020818-T01-IM6

rm(list = ls())
setwd('~/Documents/MSKCC/07_FacetsReview/')


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

TP_sample_facets = facets::readSnpMatrix(filename = 'Data4Analysis/countsMerged____P-0020818-T01-IM6_P-0020818-N01-IM6.dat.gz')
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
chr1_raw = post_process$jointseg[which(post_process$jointseg$chrom == 1), ]
chr1_raw_seg1 = chr1_raw[which(chr1_raw$seg == 1), ]
chr1_raw_hetero = chr1_raw_seg1[!is.na(chr1_raw_seg1$valor), ]

#' get the heterozygous positions:
chr1_heterozygous = TP_sample[which(TP_sample$Chromosome == 1 & TP_sample$Position %in% chr1_raw_hetero$maploc), ]
chr1_heterozygous$BAF = chr1_heterozygous$File2A / (chr1_heterozygous$File2A + chr1_heterozygous$File2R)

a = ggplot(chr1_raw_hetero, aes(x = maploc, y = valor)) + geom_point() +
        theme(aspect.ratio = 1)
b = ggplot(chr1_heterozygous, aes(x = Position, y = BAF)) + geom_point() +
        scale_y_continuous(expand = c(0,0),
                           limits = c(0,1)) +
        theme(aspect.ratio = 1)

a / b

head(chr1_heterozygous)
head(chr1_raw_hetero)

x = merge(chr1_raw_hetero, chr1_heterozygous, by.x = 'maploc', by.y = 'Position')
ggplot(x, aes(x = valor, y = BAF)) + geom_point() +
        theme(aspect.ratio = 1) +
        geom_smooth(method = "lm")


plot(chr1_heterozygous$valor)
plot(chr1_heterozygous$vafT)

cor.test(chr1_heterozygous$vafT, chr1_heterozygous$valor)

nrow(chr1_heterozygous)
head(chr1_heterozygous)
plot(chr1_raw_seg1$vafT)




plot(chr1_raw_seg1$valor)

segments(segstart, sqrt(mafR), segend, sqrt(mafR), lwd=1.75, col=2)



#' We get the following output:
facets::plotSample(x = post_process, plot.type = 'naive', sname = 'Example: P-0020818-T01-IM6')
par(mfrow = c(1,1))
facets::logRlogORspider(cncf = facets_fit$cncf, dipLogR = facets_fit$dipLogR)




# minor-allele log-odds
chrom <- out$chrom
maploc <- out$maploc
rCountN <- out$rCountN
rCountT <- out$rCountT
vafT <- out$vafT
vafN <- out$vafN
het <- out$het

#' log-odds / BAF
rcmat = round(cbind(vafT[het==1]*rCountT[het==1], 
                    (1-vafT[het==1])*rCountT[het==1], 
                    vafN[het==1]*rCountN[het==1], 
                    (1-vafN[het==1])*rCountN[het==1]))

lorvar <- valor <- rep(NA_real_, length(out$maploc))
valor[het==1] <- log(rcmat[,1]+0.5) - log(rcmat[,2]+0.5) - log(rcmat[,3]+0.5) + log(rcmat[,4]+0.5)


## illustrative example for GC-bias
plot(pctgc, log2(tcount*tscl) - log2(ncount), ylim = c(-1, 1),
     xlab = 'percent GC', ylab = 'cnlr', main = 'Example: Given cnlr within GC-bins')

lines(gcb, col = 'red', lwd = 2)
legend('bottomleft',
       col = 'red',
       lty = 1, 
       lwd = 2,
       c('Facets-default smoothing: 0.2'))


## continue with BAF ~ VALOR explanation
## extract heterozygous SNPs on X-chromosome
## look into male samples and see if there is some evidence for Y-chromosome (allelic) imbalance
## PAR1 or PAR2 region; using WES!
b$pmat
b$pmat = b$pmat[which(b$pmat$chrom == 1), ]

#' just focusing on chromosome 1; heterozygous SNPs
pmat_het = b$pmat[which(b$pmat$het == 1 & b$pmat$keep == 1), ]
#' there are 2,926 heterozygous loci on chromosome 1 which will be used for allelic copy number determination
joint_het = b$jointseg[which(b$jointseg$chrom == 1), ]
joint_het = joint_het[which(joint_het$maploc %in% pmat_het$maploc), ]

#' manually calculate the BAF in the Tumor:
#' rCountT = total coverage
#' vafT = variant-allele frequency
joint_het$BAF = joint_het$rCountT * joint_het$vafT / joint_het$rCountT
joint_het[which.max(joint_het$BAF), ]
plot(joint_het$valor, joint_het$BAF, xlab = 'logOR', ylab = 'BAF', main = 'Heterozygous SNPs on chromosome 1')
legend(legend = paste0('pearsons correlation = ', round(cor.test(joint_het$valor, joint_het$BAF)$estimate, 3)),
       x = -2.2, y = 0.88)

