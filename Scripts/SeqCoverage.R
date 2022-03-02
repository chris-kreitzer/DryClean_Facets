## Raw sequencing coverage of selected regions:
## Here, we compore Normal/Tumor samples in IGV using megadepth tool
## 
## This script should serve as aid to understand the Facets output
## 
## 03/01/2022
## chris-kreitzer

## Basically, there are two tools which are of interest to use
## https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html
## https://bioconductor.org/packages/devel/bioc/manuals/megadepth/man/megadepth.pdf

clean()
gc()
.rs.restartR()
setup(working.path = '~/Documents/MSKCC/dmp-2021/BAMs/P-0001396-T05-IM6_N01-IM6/')

# BiocManager::install("megadepth")
# megadepth::install_megadepth()
library(megadepth)


#' convert Tumor sample to BigWig file
megadepth::bam_to_bigwig(bam_file = 'DE840153-T.bam', 
                         prefix = 'DE840153-T', 
                         min_unique_qual = 30, 
                         double_count = T)

#' make sure that chromosome names are equal in BED and BigWig file
BED = data.frame(chr = c(7, 17, 'X', 10),
                 start = c(55086710, 37844347, 47004620, 89623382),
                 end = c(55279321, 37884911, 47046212, 89731687),
                 name = c('EGFR', 'ERBB2', 'RBM10', 'PTEN'))

#' make sure there are no quotes or colnames in the output .bed file
#' check the bed file structure with
#' @example bed = rtracklayer::import('4Genes.bed')
write.table(BED, file = 'genes.bed', sep = '\t', quote = F, row.names = F, col.names = F)

#' Compute the coverage of selected regions
bw_cov = get_coverage(bigwig_file = 'DE840153-T.all.bw', op = 'mean', annotation = 'genes.bed')


## bamCoverage() example
bamCoverage -b YB324274-N.bam -o Normal_PTEN.bw -of 'bigwig' -bs 1 -r 10:89623382:89731687 --minMappingQuality 30










## Read in the base-pair coverage data
if (!xfun::is_windows()) {
  regionCov <- derfinder::getRegionCoverage(
    regions = bed,
    files = example_bw,
    verbose = FALSE
  )
  ## Summarize the base-pair coverage data.
  ## Note that we have to round the mean to make them comparable.
  testthat::expect_equivalent(
    round(sapply(regionCov[c(1, 3:4, 2)], function(x) mean(x$value)), 2),
    bw_cov$score,
  )


d

a = readRDS('~/Documents/MSKCC/07_FacetsReview/DryClean/PON_BRCA/BRCA_PON_list.rds')





