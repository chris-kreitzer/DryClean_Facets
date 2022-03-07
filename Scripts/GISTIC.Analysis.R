## Prepare FACETS default runs for GISTIC 2.0 analysis ----
## Start script wrapping: 04/07/2020
## Modificatons made: 04/14/20
## currently, program run solely on R - not applicable to command line yet

prepare.GISTIC = function(path.seg.files, sample.list = NULL, output.directory, suffix, 
                          concordance = NULL, SNP.cutoff = 5, seg.mean.cutoff.deletion = -0.2,
                          seg.mean.cutoff.amplification = 0.2){
  
  # load input samples
  selected.samples = as.character(unique(sample.list)) 
  SNP.cutoff = as.numeric(SNP.cutoff) / 100
  seg.mean.cutoff.deletion = as.numeric(seg.mean.cutoff.deletion)
  seg.mean.cutoff.amplification = as.numeric(seg.mean.cutoff.amplification)
  output.directory = as.character(output.directory)
  
  print('Input SNP.cutoff is interpreted as % value; e.g. 5')
  print(paste('@least', SNP.cutoff, '% of SNP needs to be heterozygous on a segment to be included', sep = ' '))
  
  if(!is.null(path.seg.files)){
    if(!is.null(sample.list)){
      listed_files = list.files(path = path.seg.files, full.names = T) 
      
      all_files = c()
      for(samples in 1:length(selected.samples)){
        sample.input = grep(selected.samples[samples], listed_files, value = T)
        all_files = c(all_files, sample.input)
      }
    }
  }
  
  if(is.null(concordance)){
    
    ## create empty data frames
    SNP.container = data.frame()
    Segment.container = data.frame()
    
    print(length(unique(all_files)))
    
    for(files in 1:length(all_files)){
      # work on segments
      print(paste('currently @:', ' ', gsub(pattern = getwd(), replace = '', all_files[files]), sep = ''))
      seg.input = load(all_files[files])
      Gistic_segments_ICN = out$IGV
      Gistic_segments_ICN$ID = gsub(pattern = '_P.*$', replacement = '', x = Gistic_segments_ICN$ID)
      # work on SNP locations
      GISTIC.SNPs = out$jointseg[, c('chrom', 'maploc')]
      GISTIC.SNPs$Name = paste0(str_extract(all_files[files], 
                                            'P-[0-9]{7}-[A-Z]{1}[0-9]{2}-IM[0-9]{1}'), ':', 'SNP', 
                                seq(from = 1, to = nrow(GISTIC.SNPs), by = 1))
      GISTIC.SNPs = GISTIC.SNPs[,c(3,1,2)]
      colnames(GISTIC.SNPs) = c('name', 'chromosome', 'location')
      
      SNP.container = rbind(SNP.container, GISTIC.SNPs)
      Segment.container = rbind(Segment.container, Gistic_segments_ICN)
      
      rm(out)
      rm(fit)
    }
    
  } else{
    
    SNP.container = data.frame()
    Segment.container = data.frame()
    
    print(length(unique(all_files)))
    
    for(files in 1:length(all_files)){
      
      # work on segments: exclude FALSE POSITIVE
      print(paste('currently @:', ' ', gsub(pattern = getwd(), replace = '', all_files[files]), sep = ''))
      seg.input = load(all_files[files])
      Gistic_segments = out$IGV
      Gistic_segments$ID = gsub(pattern = '_P.*$', replacement = '', x = Gistic_segments$ID)
      ## prepare data from integer copy numbers; include cutoffs
      Gistic.ICN = fit$cncf
      Gistic.ICN$keep = ifelse((Gistic.ICN$nhet / Gistic.ICN$num.mark) > SNP.cutoff & !is.na(Gistic.ICN$mafR.clust), 'keep', 'remove')
      ## assign category to integer copy numbers:
      Gistic.ICN$nA = Gistic.ICN$tcn.em - Gistic.ICN$lcn.em
      Gistic.ICN$classification = ifelse(Gistic.ICN$nA == Gistic.ICN$lcn.em, 'diploid',
                                         ifelse(Gistic.ICN$nA > Gistic.ICN$lcn.em & Gistic.ICN$lcn.em != 0, 'gain', 'loss'))
      
      ## exclude FALSE positives (ICN and seg.mean)
      Gistic_segments_ICN = merge(Gistic.ICN, Gistic_segments, by.x = 'start', by.y = 'loc.start')
      
      Gistic_segments_ICN$keep2 = ifelse(Gistic_segments_ICN$classification == 'gain' & Gistic_segments_ICN$seg.mean > seg.mean.cutoff.deletion, 'keep',
                                         ifelse(Gistic_segments_ICN$classification == 'loss' & Gistic_segments_ICN$seg.mean < seg.mean.cutoff.amplification, 'keep',
                                                ifelse(Gistic_segments_ICN$classification == 'diploid' & Gistic_segments_ICN$seg.mean <= seg.mean.cutoff.amplification &
                                                         Gistic_segments_ICN$seg.mean >= seg.mean.cutoff.deletion, 'keep', 'remove')))
      
      Gistic_segments_ICN$keep2 = ifelse(is.na(Gistic_segments_ICN$keep2), 'remove', Gistic_segments_ICN$keep2)
      
      ## set segmean of false positive to 0
      Gistic_segments_ICN$seg.mean[which(Gistic_segments_ICN$keep == 'remove' | Gistic_segments_ICN$keep2 == 'remove')] = 0
      
      # subset to only those columns needed
      Gistic_segments_ICN = Gistic_segments_ICN[, c('ID', 'chrom.x', 'start', 'loc.end', 'num.mark.x', 'seg.mean')]
      colnames(Gistic_segments_ICN) = c('ID', 'chromosome', 'start', 'end', 'SNP.positions', 'seg.mean')
      
      ##################### work on SNP locations ----
      GISTIC.SNPs = out$jointseg[, c('chrom', 'maploc')]
      GISTIC.SNPs$Name = paste0(str_extract(all_files[files], 
                                            'P-[0-9]{7}-[A-Z]{1}[0-9]{2}-IM[0-9]{1}'), ':', 'SNP', 
                                seq(from = 1, to = nrow(GISTIC.SNPs), by = 1))
      GISTIC.SNPs = GISTIC.SNPs[,c(3,1,2)]
      colnames(GISTIC.SNPs) = c('name', 'chromosome', 'location')
      
      SNP.container = rbind(SNP.container, GISTIC.SNPs)
      Segment.container = rbind(Segment.container, Gistic_segments_ICN)
      
      rm(out)
      rm(fit)
      
    }
  }
  
  # prepare SNP file: remove duplicated entries
  unique.SNPs = data.frame()
  for(chromo in unique(SNP.container$chromosome)){
    print(SNP.container$chromosome[chromo])
    sub.chrom = SNP.container[SNP.container$chromosome == chromo,, drop = F]
    sub.chrom.without.duplicate = sub.chrom[!duplicated(sub.chrom$location),, drop = F]
    unique.SNPs = rbind(unique.SNPs, sub.chrom.without.duplicate)
  }
  
  dir.out = dir.create(path = paste(getwd(), '/', output.directory, sep = ''))
  
  ## output of segmentation data
  write.table(x = Segment.container, 
              file = paste(getwd(), '/', output.directory, '/', 'Segmentation.data_', as.character(suffix), '.txt', sep = ''),
              sep = '\t', row.names = F, col.names = T)
  
  ## output of SNP data
  write.table(x = unique.SNPs, 
              file = paste(getwd(), '/', output.directory, '/', 'Markerfile_SNPs_',  as.character(suffix), '.txt', sep = ''),
              sep = '\t', row.names = F, col.names = T)
  
  return(list(unique.SNPs, Segment.container))
}

## example run:
# prepare.GISTIC(path.seg.files = '~/Documents/MSKCC/ATM_project/Segmentation_Files/try2/', 
#                sample.list = ATM, output.directory = 'GISTIC.out', suffix = 'ATM_cohort')



#### INPUT GISTIC out data files and modify ----

## GISTIC output data:
## all_data_by_gene | gene level ~ patient occurrence
## del_genes_conf90 | peak regions and associated genes with q value
## amp_genes_conf90 | peak regions and associated genes with q value
## all_lesions.conf90 | all significant peak regions (DEL and AMP) ~ sample level
## 
## GISTIC file: all G scores and genomic positions

## write function which modify GISTIC output ----

GISTIC_analysis = function(path_output_files, 
                           gene.lookup = NULL, 
                           cancer.genes, 
                           suffix){
  
  
  #' load all GISTIC output files
  all_output_files = list.files(path = path_output_files, full.names = T)
  ## load cancer.genes = IMPACT pathway genes:
  load('~/Documents/MSKCC/00_Data/OncoPath12.Rdata')
  
  ## check if all_thresholded.by_genes is available
  if(all(grepl('.thresholded.by_gene.*', all_output_files) == FALSE)){
      stop(cat('all_thresholded.by_genes is missing'))
    } else {
    all_data_genes = read.csv(file = grep('.thresholded.by_gene.*', all_output_files, value = T), sep = '\t')
    }
  
  ## modify gene.lookup
  if(!is.null(gene.lookup)){
    sample_gene.lookup = all_data_genes[all_data_genes$Gene.Symbol %in% as.character(gene.lookup),, drop = F]
    sample_gene.lookup = as.data.frame(t(sample_gene.lookup))
    sample_gene.lookup = droplevels(sample_gene.lookup[-c(1,2,3),, drop = F])
    sample_gene.lookup$ID = row.names(sample_gene.lookup)
    colnames(sample_gene.lookup) = c(as.character(gene.lookup), 'ID')
    rownames(sample_gene.lookup) = NULL
    sample_gene.lookup = sample_gene.lookup[, c(2,1)]
  } 
  
  ## all_lesions.conf.XX
  if(all(grepl('.lesions.conf.*', all_output_files) == FALSE)){
    stop(cat('all_lesions.conf_XX file is missing'))
  } else {
    all_lesions = read.csv(file = grep('.lesions.conf.*', all_output_files, value = T), sep = '\t')
  }
  
  ## modify all_lesions_data
  all_lesions = all_lesions[grep('0.*', all_lesions$Amplitude.Threshold),, drop = F]
  all_lesions_subset1 = all_lesions[, c('Unique.Name', 'Descriptor', 'Peak.Limits', 'q.values')]
  all_lesions_subset2 = all_lesions[, c(10:ncol(all_lesions))]
  all_lesions_subset2$X = NULL
  
  for(i in 1:ncol(all_lesions_subset2)){
    all_lesions_subset2[, i] = ifelse(all_lesions_subset2[, i] == 2, 1,
                                      ifelse(all_lesions_subset2[, i] == 1, 1, 0))
  }
  sample.summary = data.frame(Sample.summary = rowSums(all_lesions_subset2))
  all_lesions_out = cbind(all_lesions_subset1, sample.summary)
  all_lesions_out$Peak.Limits = gsub(pattern = '.p.*', replacement = '', all_lesions_out$Peak.Limits)
  
  
  ## look into peaks: deleted and amplified genes
  if(all(grepl('.l_genes.conf.*', all_output_files) == FALSE)){
    stop(cat('file: del_genes.conf_XX is missing'))
  } else {
    del.genes = read.csv(file = grep('.l_genes.conf.*', all_output_files, value = T), sep = '\t')
  }
  
  ## modify genes in deleted region:
  IMPACT.genes = do.call('c', OncoPath_12)
  names(IMPACT.genes) = NULL
  IMPACT.genes = as.character(unique(IMPACT.genes))
  rm(OncoPath_12)
  
  del.genes = del.genes[, c(-1)]
  del.genes$X = NULL
  colnames(del.genes) = gsub('X', '', colnames(del.genes))
  del.genes = del.genes[-c(1, 2, 3),, drop = F]
  del.genes$cytoband = NULL
  
  gene.out = c()
  gene.onco = c()
  
  for(i in 1:length(del.genes)){
    gene.length.all = sum(del.genes[i] != '')
    cancer.genes = length(base::intersect(del.genes[, i], IMPACT.genes))
    gene.out = c(gene.out, gene.length.all)
    gene.onco = c(gene.onco, cancer.genes)
  }
  
  ## prepare all all genes
  chromosome = colnames(del.genes)
  del.out = data.frame(chromosome = chromosome, 
                   genes_in_peak = gene.out,
                   IMPACT_genes_peak = gene.onco)
  
  del.out$chromosome = gsub('[a-z]{1}.*', '', del.out$chromosome)
  del.out = aggregate(.~ chromosome, data = del.out, sum)
  
  missing.chr = setdiff(1:22, as.numeric(del.out$chromosome))
  missing.data = data.frame(chromosome = missing.chr,
                            genes_in_peak = rep(0, length(missing.chr)),
                            IMPACT_genes_peak = rep(0, length(missing.chr)))
  
  del.out = rbind(missing.data, del.out)
  del.out$status = rep('Deletion', nrow(del.out))
  del.out$group = rep(as.character(suffix), nrow(del.out))
  
  ## modify amplified genes
  if(all(grepl('.mp_genes.conf.*', all_output_files) == FALSE)){
    stop(cat('file: amp_genes.conf_XX is missing'))
  } else {
    amp.genes = read.csv(file = grep('.mp_genes.conf.*', all_output_files, value = T), sep = '\t')
  }
  
  ## modify amplified gene input
  amp.genes = amp.genes[, c(-1)]
  amp.genes$X = NULL
  colnames(amp.genes) = gsub('X', '', colnames(amp.genes))
  amp.genes = amp.genes[-c(1, 2, 3),, drop = F]
  amp.genes$cytoband = NULL
  
  gene.out = c()
  gene.onco = c()
  
  for(i in 1:length(amp.genes)){
    gene.length.all = sum(amp.genes[i] != '')
    cancer.genes = length(base::intersect(amp.genes[, i], IMPACT.genes))
    gene.out = c(gene.out, gene.length.all)
    gene.onco = c(gene.onco, cancer.genes)
  }
  
  ## prepare all all genes
  chromosome = colnames(amp.genes)
  amp.out = data.frame(chromosome = chromosome, 
                       genes_in_peak = gene.out,
                       IMPACT_genes_peak = gene.onco)
  
  amp.out$chromosome = gsub('[a-z]{1}.*', '', amp.out$chromosome)
  amp.out = aggregate(.~ chromosome, data = amp.out, sum)
  
  missing.chr = setdiff(1:22, as.numeric(amp.out$chromosome))
  missing.data = data.frame(chromosome = missing.chr,
                            genes_in_peak = rep(0, length(missing.chr)),
                            IMPACT_genes_peak = rep(0, length(missing.chr)))
  
  amp.out = rbind(missing.data, amp.out)
  amp.out$status = rep('Amplification', nrow(amp.out))
  amp.out$group = rep(as.character(suffix), nrow(amp.out))
  
  ## work on Gistic output
  if(all(grepl('.cores.gistic.*', all_output_files) == FALSE)){
    stop(cat('file: scores.gistic is missing'))
  } else {
    gistic = read.csv(file = grep('.cores.gistic.*', all_output_files, value = T), sep = '\t')
  }
  
  gistic$group = rep(as.character(suffix), nrow(gistic))
  
  processed.gistic = list(all_data_genes,
                          sample_gene.lookup, 
                          all_lesions_out, 
                          del.out, 
                          amp.out, 
                          gistic)
  names(processed.gistic) = c('Sample_by_gene',
                              'Altered_Samples_in_selected_gene', 
                              'significant_peak_regions', 
                              'deleted_genes_in_peak',
                              'amplified_genes_in_peak', 
                              'Gistic_estimates')
  
  return(processed.gistic)
}

#' out