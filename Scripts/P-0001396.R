## P-0001396-T05-IM6-N01
## 
## one case example: 
## Can we detect noise in the Normal and if so
## are we receiving a 'clean' tumor sample
## 
## compare raw coverage (4 selected genes)
## look into cBIO output
## run normal Facets
## run Facets with cleaned data
## 
## 03/02/2022
## chris-kreitzer

.rs.restartR()
clean()
gc()
setup(working.path = '~/Documents/MSKCC/07_FacetsReview/DryClean/')
setup(working.path = '~/Documents/GitHub/DryClean_Facets/')

## Data
PON_normalized = vroom::vroom('DataProcessed/PON_normalized.txt')
PON_normalized = as.data.frame(PON_normalized)
row.names(PON_normalized) = PON_normalized$duplication
PON_normalized$duplication = NULL


## Population-wise: which bins show the greatest variability:
average = apply(PON_normalized, 1, FUN = mean)
variance = apply(PON_normalized, 1, FUN = var)

#' which bins are overrepresented? cluster around which chromosomes, genes and locations?
#' here, we are only asking for >99% percentile
#' write a little function that fetches regions that are over represented in the NORMALS

Seq_overRepresentation = function(sampleID, 
                                  selection_criterion = 'average', 
                                  quantile_selected = 0.99,
                                  break_threshold = 1000000){
  message('This function currently just works with PON_normalized data structure')
  message('Currently, 1Mb or 1,000 000 bps are required to find a new block')
  
  if(selection_criterion == 'average'){
    message('apply(PON_normalized, 1, FUN = mean) was run on population; n = 1,000 samples')
  }
  
  #' select bins with >> average (1,000) seq. depth
  #stopifnot(exprs = 'average' %in% ls(environment()))
  bin_overrepresented = names(sort(average[which(average > quantile(x = average, probs = quantile_selected))], decreasing = T))
  data_sub = PON_normalized[, grepl(pattern = sampleID, x = colnames(PON_normalized))]
  data_sub$bin = row.names(data_sub)
  row.names(data_sub) = NULL
  data_sub = data_sub[which(data_sub$bin %in% bin_overrepresented), ]
  print(grep(pattern = sampleID, x = colnames(PON_normalized), value = T))
  

  #' split bin and order data frame
  data_sub$Chromosome = unlist(strsplit(as.character(data_sub$bin), ';'))[2*(1:length(data_sub$bin)) - 1]
  data_sub$Position = unlist(strsplit(as.character(data_sub$bin), ';'))[2*(1:length(data_sub$bin))]
  data_sub = data_sub[order(data_sub[, "Chromosome"], data_sub[, "Position"]), ]
  
  data_sub$diff = NA
  data_sub$Position = as.numeric(as.character(data_sub$Position))
  
  #' insert the breakpoints
  i = 1
  while(i < nrow(data_sub)){
    data_sub$diff[i] = data_sub$Position[i+1] - data_sub$Position[i]
    i = i + 1
  }
  
  data_sub$breaks = NA
  data_sub$breaks = ifelse(abs(data_sub$diff) > break_threshold, 'end', NA)
  
  #' find the start positions
  ii_end = which(data_sub$breaks == 'end')
  ii_start = ii_end + 1
  data_sub$breaks[ii_start] = 'start'
  
  #' exclude NAs
  data_sub$breaks[1] = 'start'
  data_sub = data_sub[!is.na(data_sub$breaks), ]
  
  return(data_sub)
  
}

y = Seq_overRepresentation(sampleID = 'P-0001396-T05-IM6.*', selection_criterion = 'average')

#' create .bed file of overrepresented regions:
Normal_over_bed = data.frame()
for(i in seq(1, nrow(y), 2)){
  chr = y$Chromosome[i]
  start = y$Position[i]
  end = y$Position[i+1]
  sub = data.frame(chromosome = chr,
                   start = start,
                   end = end)
  Normal_over_bed = rbind(Normal_over_bed, sub)
}
Normal_over_bed = Normal_over_bed[-nrow(Normal_over_bed), ]
Normal_over_bed$chromosome = as.numeric(as.character(Normal_over_bed$chromosome))
Normal_over_bed = Normal_over_bed[order(Normal_over_bed[, "chromosome"], Normal_over_bed[, "start"]), ]
ii = which(Normal_over_bed$end < Normal_over_bed$start)
Normal_over_bed = Normal_over_bed[-ii, ]

write.table(Normal_over_bed, file = 'DataProcessed/normal.bed', 
            sep = '\t', quote = F, col.names = F, row.names = F)


#' human hg19_genes:
human_hg19_genes = vroom::vroom('~/Documents/MSKCC/dmp-2021/human_hg19_genes.tsv.gz')
human_hg19_genes = as.data.frame(human_hg19_genes)
human_hg19_genes = human_hg19_genes[, c('name', 'chrom', 'txStart', 'txEnd', 'name2')]

human_genes = data.frame()
for(i in unique(human_hg19_genes$name2)){
  print(i)
  if(nrow(human_hg19_genes[which(human_hg19_genes$name2 == i), ]) > 1){
    data_select = human_hg19_genes[which(human_hg19_genes$name2 == i), ]
    data_select = data_select[1, ]
  } else {
    data_select = human_hg19_genes[which(human_hg19_genes$name2 == i), ]
  }
  human_genes = rbind(human_genes, data_select)
}


human_genes.bed = human_genes[, c('chrom', 'txStart', 'txEnd', 'name2')]
human_genes.bed$chrom = sub(pattern = '^chr', replacement = '', human_genes.bed$chrom)
human_genes.bed$chrom = as.numeric(as.character(human_genes.bed$chrom))
human_genes.bed = human_genes.bed[order(human_genes.bed[, "chrom"], human_genes.bed[, "txStart"]), ]

#' just IMPACT genes
human_genes.bed = human_genes.bed[which(human_genes.bed$name2 %in% IMPACT468), ]
write.table(human_genes.bed, file = 'DataProcessed/hg19_genes.bed', quote = F, sep = '\t', row.names = F, col.names = F)

##-----------------------------------------------------------------------------
## look into the GRanges object of the NORMAL P-0001396-T05-IM6/N01
normal_table = readRDS('PON_BRCA/normal_table.rds')
tumor_table = readRDS('TUMOR_BRCA/tumor_table.rds')
GrNormal = readRDS('~/Documents/MSKCC/07_FacetsReview/DryClean/PON_BRCA/sample771.rds')
GrTumor = readRDS('~/Documents/MSKCC/07_FacetsReview/DryClean/TUMOR_BRCA/sample3.rds')

#' modify Normal
Normal_raw = as.data.frame(GrNormal)
Normal_raw$bin = paste(Normal_raw$seqnames, Normal_raw$start, sep = ';')
Normal_raw$end = NULL
Normal_raw$strand = NULL
Normal_raw$width = NULL
Normal_raw$seqnames = as.character(as.factor(Normal_raw$seqnames))
Normal_raw$seqnames[which(Normal_raw$seqnames == 'X')] = 23
Normal_raw$seqnames[which(Normal_raw$seqnames == 'Y')] = 24
Normal_raw$seqnames = as.numeric(as.character(Normal_raw$seqnames))
Normal_raw$start = as.numeric(as.integer(Normal_raw$start))
Normal_raw = Normal_raw[order(Normal_raw[, "seqnames"], Normal_raw[, "start"]), ]

#' modify Tumor
Tumor_raw = as.data.frame(GrTumor)
Tumor_raw$bin = paste(Tumor_raw$seqnames, Tumor_raw$start, sep = ';')
Tumor_raw$end = NULL
Tumor_raw$strand = NULL
Tumor_raw$width = NULL
Tumor_raw$seqnames = as.character(as.factor(Tumor_raw$seqnames))
Tumor_raw$seqnames[which(Tumor_raw$seqnames == 'X')] = 23
Tumor_raw$seqnames[which(Tumor_raw$seqnames == 'Y')] = 24
Tumor_raw$seqnames = as.numeric(as.character(Tumor_raw$seqnames))
Tumor_raw$start = as.numeric(as.integer(Tumor_raw$start))
Tumor_raw = Tumor_raw[order(Tumor_raw[, "seqnames"], Tumor_raw[, "start"]), ]


#' Look at specific genes:
#' PTEN
PTEN = exonic_structure(gene = 'PTEN', type = 'exons')
PTEN = PTEN[which(PTEN$transcript_id == 'ENST00000371953.3'), ]

start = min(PTEN$start)
end = max(PTEN$end)

Norm_pten = Normal_raw[which(Normal_raw$seqnames == 10 & Normal_raw$start >= start & Normal_raw$start <= end), ]
Tumor_pten = Tumor_raw[which(Tumor_raw$seqnames == 10 & Tumor_raw$start >= start & Tumor_raw$start <= end), ]

#' Visualization:
pten_raw_merged = merge(Norm_pten, Tumor_pten, by = 'bin')
ggplot(pten_raw_merged, aes(x = reads.corrected.x, y = reads.corrected.y)) +
  geom_density2d_filled(alpha = 0.5) +
  geom_point(size = 0.3) +
  scale_color_viridis_c() +
  theme(aspect.ratio = 1) +
  geom_abline(intercept = 0, slope = 1)


par(mfrow = c(2,1))
plot(Norm_pten$reads.corrected, xaxt = 'n', xlab = '')
plot(Tumor_pten$reads.corrected, xaxt = 'n', xlab = '')

head(Normal_raw)
head(Tumor_raw)
str(Normal_raw)
levels(Tumor_raw$seqnames)
head(Normal_raw)

#' Look into EGFR, where we do know that we have a gross deviation
#' PTEN
EGFR = exonic_structure(gene = 'EGFR', type = 'exons')
EGFR = EGFR[which(EGFR$transcript_id == 'ENST00000455089.1'), ]

start = min(EGFR$start)
end = max(EGFR$end)
chrom = 7

Norm_pten = Normal_raw[which(Normal_raw$seqnames == chrom & Normal_raw$start >= start & Normal_raw$start <= end), ]
Tumor_pten = Tumor_raw[which(Tumor_raw$seqnames == chrom & Tumor_raw$start >= start & Tumor_raw$start <= end), ]

#' Visualization:
pten_raw_merged = merge(Norm_pten, Tumor_pten, by = 'bin')
ggplot(pten_raw_merged, aes(x = reads.corrected.x, y = reads.corrected.y)) +
  geom_density2d_filled(alpha = 0.5) +
  geom_point(size = 0.3) +
  scale_color_viridis_c() +
  theme(aspect.ratio = 1) +
  geom_abline(intercept = 0, slope = 1)


par(mfrow = c(2,1))
plot(Norm_pten$reads.corrected, xaxt = 'n', xlab = '')
plot(Tumor_pten$reads.corrected, xaxt = 'n', xlab = '')








a = FacetsDC::run_facets_cleaned(read_counts = '~/Desktop/countsMerged____P-0001396-T05-IM6_P-0001396-N01-IM6.dat.gz',
                                 read_cleaned = '~/Desktop/sample4.rds.rds', MODE = 'union')




x1 = facets::readSnpMatrix(filename = '~/Desktop/countsMerged____P-0001396-T05-IM6_P-0001396-N01-IM6.dat.gz')
x2 = facets::preProcSample(rcmat = x1)
x3 = facets::procSample(x2, cval = 150)
x4 = facets::emcncf(x3)
x4$loglik

length(a$segs$chrom[which(a$segs$tcn.em == 2 & a$segs$lcn.em != 1)])
View(a$segs)

cbio = read.csv('~/Downloads/mskimpact_segments (1).seg', sep = '\t')

##-----------------------------------------------------------------------------
paths = read.csv('~/Documents/MSKCC/07_FacetsReview/DryClean/DataProcessed/PON_BRCA_Paths.txt', sep = '\t', header = F)
View(paths)
normalized = data.table::fread(input = '~/Documents/MSKCC/07_FacetsReview/DryClean/Data4Analysis/normalizedPON.txt', sep = '\t')
sample_norm = normalized[, grep(pattern = '^P-0001396-T05-IM6.*', x = colnames(normalized), value = T)]
head(sample_norm)
which(sample_norm) == T
head(colnames(normalized))
grep(pattern = '^P-0001396.*', x = colnames(normalized), value = T)



if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("megadepth")


setwd('~/Documents/MSKCC/dmp-2021/BAMs/P-0001396-T05-IM6_N01-IM6/')
megadepth::bam_to_bigwig(bam_file = 'DE840153-T.bam', prefix = 'TEST.bg', min_unique_qual = 30, double_count = T)





old = vroom::vroom('~/Documents/MSKCC/07_FacetsReview/DryClean/Data4Analysis/normalizedPON.txt')
new = vroom::vroom('~/Documents/MSKCC/07_FacetsReview/DryClean/DataProcessed/PON_normalized.txt')


View(old[1:5, 1:5])
View(new[1:5, 1:5])

u = old[,c(1:3)]
r = new[,c(1:3)]

View(u)
View(r)


n = merge(u, r, by = 'duplication')
View(n)





