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
chrom = 10

EGFR = exonic_structure(gene = 'EGFR', type = 'exons')
EGFR = EGFR[which(EGFR$transcript_id == 'ENST00000455089.1'), ]
start = min(EGFR$start)
end = max(EGFR$end)
chrom = 7

ERBB2 = exonic_structure(gene = 'ERBB2', type = 'exons')
ERBB2 = ERBB2[which(ERBB2$transcript_id == 'ENST00000584601.1'), ]
start = min(ERBB2$start)
end = max(ERBB2$end)
chrom = 17

RBM10 = exonic_structure(gene = 'RBM10', type = 'exons')
RBM10 = RBM10[which(RBM10$transcript_id == 'ENST00000377604.3'), ]
start = min(RBM10$start)
end = max(RBM10$end)
chrom = 23

#' function
raw_gene_seq_vis = function(normal, tumor, gene, chrom, start, end){
  Norm = normal[which(normal$seqnames == chrom & normal$start >= start & normal$start <= end), ]
  Tumor = tumor[which(tumor$seqnames == chrom & tumor$start >= start & tumor$start <= end), ]
  
  data_merge = merge(Norm, Tumor, by = 'bin')
  
  #' Visualization
  plot = ggplot(data_merge, aes(x = reads.corrected.x, y = reads.corrected.y)) +
    stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, max(data_merge$reads.corrected.x))) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, max(data_merge$reads.corrected.y))) +
    scale_fill_distiller(palette = 'Greens', direction = 1) +
    geom_abline(intercept = 0, slope = 1) +
    theme(aspect.ratio = 1,
          panel.border = element_rect(fill = NA, size = 1.2),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = 'none',
          legend.key.size = unit(1.0, "cm"),
          legend.key.width = unit(0.8,"cm"),
          legend.key.height = unit(0.35,"cm"),
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
    labs(x = 'Normal [Seq. Coverage]', y = 'Tumor [Seq. Coverage]', title = gene)
  
  return(plot)
  
}

PTEN = raw_gene_seq_vis(normal = Normal_raw, tumor = Tumor_raw, gene = 'PTEN', chrom = 10, start = start, end = end)
EGFR = raw_gene_seq_vis(normal = Normal_raw, tumor = Tumor_raw, gene = 'EGFR', chrom = chrom, start = start, end = end)
ERBB2 = raw_gene_seq_vis(normal = Normal_raw, tumor = Tumor_raw, gene = 'ERBB2', chrom = chrom, start = start, end = end)
RBM10 = raw_gene_seq_vis(normal = Normal_raw, tumor = Tumor_raw, gene = 'RBM10', chrom = chrom, start = start, end = end)

library(cowplot)
plot_grid(PTEN, EGFR, ERBB2, RBM10, nrow = 2, ncol = 2)
## export as 7 x 7 (portrait)


##-----------------------------------------------------------------------------
## We see noise in the normal sample: Especially in those regions where we have
## more than average seq. coverage (look into the normal bed); we have 96 regions
## with aberrant signals; do we see whether those are filtered out in the decomposed 
## tumor? Inspect some regions in the normal:

## fetch regions of the normal, where we see a lot of noise;
## noise in terms of excessive seq. coverage; we have extracted those regions already:
## Now, we nee to annotate those - and see which regions are 'unspecific'; unspecific in terms of
## length, number of genes and number of IMPACT468 genes

genecode = rtracklayer::import('~/Documents/GitHub/DryClean_Facets/gencode.v19.chr_patch_hapl_scaff.annotation.gtf.gz')

.annotate_noise = function(gencode = genecode, BED.file){
  genecode = as.data.frame(genecode)
  file_2_annotate = read.csv(BED.file, sep = '\t', header = F)
  ii = which(is.na(file_2_annotate$V1))
  file_2_annotate$V1 = paste0('chr', file_2_annotate$V1)
  file_2_annotate$V1[ii] = 'chrX'
  
  annotation_out = data.frame()
  for(i in 1:nrow(file_2_annotate)){
    chrom = file_2_annotate$V1[i]
    start = file_2_annotate$V2[i]
    end = file_2_annotate$V3[i]
    length = end - start
    data_sub = genecode[which(genecode$seqnames == chrom & genecode$start >= start & genecode$end <= end), ]
    data_sub = data_sub[!duplicated(data_sub$gene_name), ]
    IMPACT_intersect = length(intersect(data_sub$gene_name, IMPACT468))
    IMPACT_genes = intersect(IMPACT468, data_sub$gene_name)
    IMPACT_genes = as.character(paste(IMPACT_genes, sep = ';', collapse = ';'))
    
    if(length(IMPACT_intersect) == 0){
      IMPACT_intersect = 0
      IMPACT_genes = NA
    }
    
    data_out = data.frame(chrom = chrom,
                          start = start,
                          end = end,
                          length = length,
                          IMPACT_intersect = IMPACT_intersect,
                          IMPACT_genes = IMPACT_genes)
    annotation_out = rbind(annotation_out, data_out)
  }
  return(annotation_out)
}

u = .annotate_noise(BED.file = '~/Documents/MSKCC/07_FacetsReview/DryClean/DataProcessed/normal.bed')
# write.table(x = u, file = 'DataProcessed/normal_bed_annotated.txt', sep = '\t', quote = F)


##-----------------------------------------------------------------------------
## compare raw tumor and cleaned tumor at selected loci from the annotated file above
cov = readRDS('~/Documents/MSKCC/07_FacetsReview/DryClean/TUMOR_BRCA/sample3.rds')
sample_clean = dryclean::start_wash_cycle(cov = cov, 
                                          mc.cores = 1, 
                                          detergent.pon.path = '~/Documents/GitHub/DryClean_Facets/detergent.rds')
Tumor_clean = as.data.frame(sample_clean)
Tumor_clean$bin = paste(Tumor_clean$seqnames, Tumor_clean$start, sep = ';')
colnames(Tumor_clean)[8] = 'reads.corrected'


#' merge the data:
data_merged = merge(Normal_raw, Tumor_clean, by = 'bin', all.y = T)

#' look into specific regions:
#' chromosome17; centromeric region: ~7 Mbp in length
chrom = 17
start = 21544504
end = 29095900

region1 = data_merged[which(data_merged$seqnames.x == chrom & 
                              data_merged$start.x >= start & 
                              data_merged$start.x <= end), ]


ggplot(region1, aes(x = reads.corrected.x, y = background)) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_distiller(palette = 'Greens', direction = 1) +
  theme(aspect.ratio = 1,
        panel.border = element_rect(fill = NA, size = 1.2),
        axis.ticks = element_blank(),
        legend.position = 'none',
        legend.key.size = unit(1.0, "cm"),
        legend.key.width = unit(0.8,"cm"),
        legend.key.height = unit(0.35,"cm"),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
  


cor.test(region1$reads.corrected.x, region1$background)


#' chromosome6; 6q arm; 11.8 Mb in length
chrom = 6
start = 117661000
end = 129478900

region2 = data_merged[which(data_merged$seqnames.x == chrom & 
                              data_merged$start.x >= start & 
                              data_merged$start.x <= end), ]


ggplot(region2, aes(x = reads.corrected.x, y = background)) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 2.2)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2.2)) +
  geom_abline(intercept = 0, slope = 1) +
  scale_fill_distiller(palette = 'Greens', direction = 1) +
  annotate(geom = 'text', x = 0.3, y = 2, 
           label = paste0('Correlation:\n ', 
                          round(cor.test(region2$reads.corrected.x, region2$background)[[4]], 3))) +
  theme(aspect.ratio = 1,
        panel.border = element_rect(fill = NA, size = 1.2),
        axis.ticks = element_blank(),
        legend.position = 'none',
        legend.key.size = unit(1.0, "cm"),
        legend.key.width = unit(0.8,"cm"),
        legend.key.height = unit(0.35,"cm"),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
  labs(x = 'normalized coverage - NORMAL', y = 'DryClean background', title = 'Chrom 6q ~ 11.8Mb')




##-----------------------------------------------------------------------------
## GISTIC output analysis:
source(file = '~/Documents/GitHub/DryClean_Facets/Scripts/GISTIC.Analysis.R')
x = GISTIC_analysis(path_output_files = '~/Documents/MSKCC/07_FacetsReview/DryClean/GISTIC_out/416318/',
                    gene.lookup = 'PTEN', suffix = 'DryClean')

y = x$Sample_by_gene[,c(1, which(colnames(x$Sample_by_gene) == 'P.0001396.T05.IM6'))]

#' here we can inspect the sample on a gene level if neccessary









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





