##-----------------------------------------------------------------------------
## Focality of copy number calls: Focus on ERBB2
## ----------------------------------------------------------------------------
##
## start: 04/27/2022
## chris-kreitzer


clean()
gc()

source('Scripts/CnLR_plot.R')

MasterFile = readRDS('~/Documents/MSKCC/07_FacetsReview/DryClean/TUMOR_BRCA/tumor_table.rds')
countfiles = list.files('~/Documents/MSKCC/07_FacetsReview/DryClean/Tumor_Countfiles/', full.names = T)
cleanedfiles = list.files('~/Documents/MSKCC/07_FacetsReview/DryClean/TUMOR_DECOMPOSED/', full.names = T)


ERBB2 = data.frame()
gene.start = 37842337
gene.end = 37886915
cval = 50
template = matrix(data = NA, nrow = 1, ncol = 17)
colnames(template) = c('chrom','seg','num.mark','nhet','cnlr.median',
'mafR','segclust','cnlr.median.clust','mafR.clust','Start','End','cf.em',
'tcn.em','lcn.em','cf','tcn','lcn')
template = as.data.frame(template)


for(i in 1:nrow(MasterFile)){
  try({
    print(MasterFile$original[i])
    count = facetsSuite::read_snp_matrix(input_file = grep(pattern = MasterFile$original[i], x = countfiles, value = T))
    count = count[, c(1,2,3,5,4,6)]
    facets_out = facetsSuite::run_facets(read_counts = count, 
                                         cval = cval, 
                                         genome = 'hg19', 
                                         seed = 100)
    segs_original = facets_out$segs[which(facets_out$segs$chrom == 17), ]
    colnames(segs_original)[10] = 'Start'
    colnames(segs_original)[11] = 'End'
    
    ERBB2_segs = segs_original[which(segs_original$Start <= gene.start & segs_original$End <= gene.end & segs_original$End >= gene.start |
                                       segs_original$Start <= gene.end & segs_original$Start >= gene.start & segs_original$End >= gene.end |
                                       gene.start <= segs_original$Start & gene.end >= segs_original$End |
                                       segs_original$Start <= gene.start & segs_original$End >= gene.end), ]
      
    fit_original = facets_fit_qc(facets_output = facets_out)$facets_qc
    
    if(nrow(ERBB2_segs) != 0){
      ERBB2_segs$sample = MasterFile$original[i]
      ERBB2_segs$qc = fit_original
      ERBB2_segs$algorithm = 'facets_original'
      
    } else {
      ERBB2_segs = template
      ERBB2_segs$sample = MasterFile$original[i]
      ERBB2_segs$qc = fit_original
      ERBB2_segs$algorithm = 'facets_original'
    }

    ##-------------------------------------------
    ## modified version
    cleaned = FacetsDC::run_facets_cleaned(read_counts = grep(pattern = MasterFile$original[i], x = countfiles, value = T),
                                           read_cleaned = grep(pattern = paste0(MasterFile$sample[i], '.rds.*'), x = cleanedfiles, value = T),
                                           MODE = 'union', 
                                           cval = cval, 
                                           seed = 100)
    
    segs_cleaned = cleaned$segs[which(cleaned$segs$chrom == 17), ]
    colnames(segs_cleaned)[10] = 'Start'
    colnames(segs_cleaned)[11] = 'End'
    
    ERBB2_segs_cleaned = segs_cleaned[which(segs_cleaned$Start <= gene.start & segs_cleaned$End <= gene.end & segs_cleaned$End >= gene.start |
                                               segs_cleaned$Start <= gene.end & segs_cleaned$Start >= gene.start & segs_cleaned$End >= gene.end |
                                               gene.start <= segs_cleaned$Start & gene.end >= segs_cleaned$End |
                                               segs_cleaned$Start <= gene.start & segs_cleaned$End >= gene.end), ]
    
    fit_cleaned = facets_fit_qc(facets_output = cleaned)$facets_qc
    
    if(nrow(ERBB2_segs_cleaned) != 0){
      ERBB2_segs_cleaned$sample = MasterFile$original[i]
      ERBB2_segs_cleaned$qc = fit_cleaned
      ERBB2_segs_cleaned$algorithm = 'dryclean'
      
    } else {
      ERBB2_segs_cleaned = template
      ERBB2_segs_cleaned$sample = MasterFile$original[i]
      ERBB2_segs_cleaned$qc = fit_cleaned
      ERBB2_segs_cleaned$algorithm = 'dryclean'
    }
    
    rm(fit_cleaned, cleaned, fit_original, count, facets_out, segs_original, segs_cleaned)
    
    ERBB2 = rbind(ERBB2, ERBB2_segs, ERBB2_segs_cleaned)
  })
  
}

write.table(ERBB2, file = 'Data_out/ERBB2_focality.txt', sep = '\t', row.names = F, quote = F)


##-----------------------------------------------------------------------------
## Assessment of ERBB2:
##-----------------------------------------------------------------------------
qc = ERBB2[which(ERBB2$qc == T), ]
qc = qc[which(qc$tcn.em != 2 & qc$lcn.em != 1), ]

# 'P-0034349-T01-IM6'
count = facetsSuite::read_snp_matrix(input_file = grep(pattern = 'P-0034349-T01-IM6', x = countfiles, value = T))
count = count[, c(1,2,3,5,4,6)]
facets_out = facetsSuite::run_facets(read_counts = count, 
                                     cval = cval, 
                                     genome = 'hg19', 
                                     seed = 100)
facets_fit_qc(facets_output = facets_out)
a = cnlr_plot(facets_out, genome = hg19)



cleaned = FacetsDC::run_facets_cleaned(read_counts = grep(pattern = 'P-0034349-T01-IM6.*', x = countfiles, value = T),
                                       read_cleaned = grep(pattern = 'sample668.rds.*', x = cleanedfiles, value = T),
                                       MODE = 'union', 
                                       cval = cval, 
                                       seed = 100)

facets_fit_qc(cleaned)
b = cnlr_plot(cleaned, genome = hg19)
b = b + labs(y = 'cleaned CnLR')

library(patchwork)
a/b


## Example 2:
# P-0027273-T01-IM6
count = facetsSuite::read_snp_matrix(input_file = grep(pattern = 'P-0027273-T01-IM6', x = countfiles, value = T))
count = count[, c(1,2,3,5,4,6)]
facets_out = facetsSuite::run_facets(read_counts = count, 
                                     cval = cval, 
                                     genome = 'hg19', 
                                     seed = 100)
facets_fit_qc(facets_output = facets_out)
a = cnlr_plot(facets_out, genome = hg19)



cleaned = FacetsDC::run_facets_cleaned(read_counts = grep(pattern = 'P-0027273-T01-IM6.*', x = countfiles, value = T),
                                       read_cleaned = grep(pattern = 'sample531.rds.*', x = cleanedfiles, value = T),
                                       MODE = 'union', 
                                       cval = cval, 
                                       seed = 100)

facets_fit_qc(cleaned)
b = cnlr_plot(cleaned, genome = hg19)
b = b + labs(y = 'cleaned CnLR')

library(patchwork)
a/b
