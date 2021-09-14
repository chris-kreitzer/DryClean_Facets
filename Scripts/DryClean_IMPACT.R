## Running the decompostion on tumor samples
## 


library(GenomicRanges)
#' read the rds object - than convert to GRanges object - then run it
#' fetch one converted sample from juno script;

cov = readRDS('~/Desktop/countsMerged____P-0002273-T01-IM3_P-0002273-N01-IM3.dat.gz.rds')

#' convert object
a = GRanges(as.data.frame(GRangesList(cov)))

cov_out = dryclean::start_wash_cycle(cov = a,
                                     detergent.pon.path = '~/Documents/MSKCC/07_FacetsReview/PON_BRCA/detergent.rds',
                                     whole_genome = F,
                                     chr = NA,
                                     germline.filter = FALSE, mc.cores = 4)



# Loading PON a.k.a detergent from path provided
# Let's begin, this is whole exome/genome
# Initializing wash cycle
# Using the detergent provided to start washing
# lambdas calculated
# calculating A and B
# calculating v and s
# Error in m.vec - s : non-conformable arrays



cov
cov_out = start_wash_cycle(cov = coverage_file, 
                           detergent.pon.path = "~/git/dryclean/inst/extdata", 
                           whole_genome = TRUE, 
                           chr = NA, 
                           germline.filter = TRUE, 
                           germline.file = "~/git/dryclean/inst/extdata/germline.markers.rds")


