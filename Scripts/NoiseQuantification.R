#' Noise detection; population wise (selected cohort)
#' Discriminate between high and low purity samples
#' 
#' n = 65 selected samples
#' 
#' 10/11/2021: adapt ideas from single run and apply to all samples
#' chris-kreitzer
#' 

rm(list = ls())
.rs.restartR()
set.seed(99)
setwd('~/Documents/GitHub/DryClean_Facets/')
source('Scripts/UtilityFunctions.R')

## Libraries and Input
library(DNAcopy)
library(facets)
library(data.table)
library(ggplot2)
library(patchwork)
library(cowplot)

DataIn = readRDS('Data_out/BRCA_workingCohort_MSK.rds')

#' probes to investigate
RecognizedProbes = as.data.frame(readRDS('Normal_samples/countsMerged____P-0000584-T03-IM6_P-0000584-N01-IM6.dat.gz.rds'))
colnames(RecognizedProbes)[1] = 'chromosome'
colnames(RecognizedProbes)[2] = 'start'
RecognizedProbes$probes = paste(RecognizedProbes$chromosome, RecognizedProbes$start, sep = ';')

#' ROI
ROI = data.frame(chromosome = c(17, 8, 1),
                 start = c(27229517, 1300569, 131354907),
                 end = c(48441764, 44119298, 239198818),
                 gene = c('ERBB2', '8p_arm', '1q_arm'))



#' TN metrics; naive log(Tumor-depth / Normal-depth)
#' @param: data = source path
#' @details: files = list.files('Tumor_countsFile/', full.names = T, all.files = F, include.dirs = T) 

TN_metrics = function(data){
  file_name = basename(data)
  countdata = facets::readSnpMatrix(data, err.thresh = 10, del.thresh = 10)
  countdata$probes = paste(countdata$Chromosome, countdata$Position, sep = ';')
  countdata = countdata[which(countdata$Chromosome == ROI[1, 'chromosome']), ]
  countdata = countdata[which(countdata$probes %in% RecognizedProbes$probes[which(RecognizedProbes$chromosome == ROI[1, 'chromosome'])]), ]
  countdata = countdata[which(countdata$Position > ROI[1, 'start'] &
                                countdata$Position < ROI[1, 'end']), ]
  countdata$TN_ratio = log(countdata$TUM.DP / countdata$NOR.DP)
  countdata$indx = seq(from = 1, to = nrow(countdata), by = 1)
  countdata$sample = file_name
  countdata$dispersion = dlrs(x = countdata$TN_ratio)
  countdata
}

# x = lapply(files, function(x) TN_metrics(x))
# y = data.table::rbindlist(x)


#' Segmentation; CBS on T/N samples
#' @param: data = source path

#' It implements our methodology for finding change-points in these data (Olshenet al., 2004), 
#' which are points after which the (log) test over reference ratios have changed location. 
#' Our model is that the change-points correspond to positions where the underlying 
#' DNA copy number has changed.

CBS_segmentation = function(data){
  name = unique(as.character(data$sample))
  CBS_data = data[,c('Chromosome', 'Position', 'TN_ratio', 'sample')]
  CNA.object = CNA(genomdat = cbind(CBS_data$TN_ratio),
                   chrom = CBS_data$Chromosome, 
                   maploc = CBS_data$Position, 
                   data.type = 'logratio', 
                   sampleid = name)
  smoothed_object = smooth.CNA(x = CNA.object)
  segmented_object = DNAcopy::segment(x = smoothed_object, 
                                      min.width = 5,
                                      undo.splits = 'sdundo', 
                                      undo.SD = 3)
  segmented_object
}


#' @import: xjuno output (from server)
TN_raw = read.csv('Data_out/TN_raw_out.txt', sep = '\t')
TN_cbs = readRDS('Data_out/TN_CBS.rds')
ERBB2_coordinates = TN_raw$indx[which(TN_raw$sample == 'countsMerged____P-0000584-T03-IM6_P-0000584-N01-IM6.dat.gz' &
                                   TN_raw$Position >= 37842338 &
                                   TN_raw$Position <= 37886915)]

#' Visualisation of TN samples and noise quantifications
plot_samples = function(data_raw, data_cbs){
  #' raw log T/N distribution
  TN_raw = ggplot(data_raw, aes(x = indx, y = TN_ratio)) +
    geom_point() +
    theme(aspect.ratio = 0.5,
          axis.line.y = element_line(colour = 'black', size = 0.2),
          panel.background = element_blank()) +
    scale_y_continuous(limits = c(-2, 2)) +
    geom_hline(yintercept = seq(-2, 2, 1), size = 0.1, linetype = 'dashed') +
    labs(title = paste0('chromosome_17q', '; dlrs: ', round(data_raw$dispersion, 3))) +
    theme(plot.margin = margin(0, 0, 0, 0, "cm")) +
    geom_segment(aes(x = ERBB2_coordinates[1], xend = ERBB2_coordinates[length(ERBB2_coordinates)],
                     y = 1.3, yend = 1.3), color = 'red', size = 1.2) +
    geom_text(x = ERBB2_coordinates[1], y = 1.45, label = "ERBB2", vjust = 'middle')
  
  #' make the histogram
  TN_dispersion = ggplot(data_raw, aes(x = TN_ratio)) + 
    geom_histogram(aes(y = ..density..), bins = 50, colour="black", fill="white") +
    labs(title = paste0('median = ', round(median(data_raw$TN_ratio), 3), '; sd = ', 
                        round(sd(data_raw$TN_ratio), 3)),
         x = 'log T/N ratio')
  
  
  #' CBS visualization:
  cbs_plot_raw = data_cbs$data
  name = substr(x = names(cbs_plot_raw)[3], start = 17, stop = 33)
  
  genomdat = cbs_plot_raw[[3]]
  maploc = 1:length(genomdat)
  segres = data_cbs$output
  
  p2 = function(){
    par(
      mar = c(4, 2, 4, 2),
      mgp = c(2, 1, 0)
    )
    plot(maploc, genomdat, 
         col = 'black', 
         pch = '.', main = name,
         ylim = c(-2, 2),
         cex = 2)
    abline(h = seq(-2, 2, 1), lty = 'dashed', lwd = 0.2)
    ii = cumsum(c(0, segres$num.mark))
    mm = segres$seg.mean
    kk = length(ii)
    segments(maploc[ii[-kk] + 1], segres$seg.mean, 
             x1 = maploc[ii[-1]], y1 = segres$seg.mean, 
             col = 'red')
  }
  
  #' make the output
  p1 = TN_raw + TN_dispersion
  plot_grid(p1, ggdraw(p2))
  
}

#' array plot for data frame()
plot_list_tn = list()
for(i in unique(TN_raw$sample)){
  data_raw = TN_raw[which(TN_raw$sample == i), ]
  data_cbs = TN_cbs[[which(names(TN_cbs) == i)]]
  plot_list_tn[[i]] = plot_samples(data_raw = data_raw, data_cbs = data_cbs)
}



###############################################################################
#' second panel: concentrate on Facets output
Facets_metrics = function(data){
  output = list()
  file_name = basename(data)
  countdata = facets::readSnpMatrix(data, err.thresh = 10, del.thresh = 10)
  countdata = facets::preProcSample(rcmat = countdata, snp.nbhd = 250)
  countdata_raw = countdata$jointseg
  countdata_raw$probes = paste(countdata_raw$chrom, countdata_raw$maploc, sep = ';')
  countdata_raw = countdata_raw[which(countdata_raw$chrom == ROI[1, 'chromosome']), ]
  countdata_raw = countdata_raw[which(countdata_raw$probes %in% RecognizedProbes$probes[which(RecognizedProbes$chromosome == ROI[1, 'chromosome'])]), ]
  countdata_raw = countdata_raw[which(countdata_raw$maploc > ROI[1, 'start'] &
                                        countdata_raw$maploc < ROI[1, 'end']), ]
  countdata_raw$TN_ratio = countdata_raw$cnlr
  countdata_raw$indx = seq(from = 1, to = nrow(countdata_raw), by = 1)
  countdata_raw$sample = file_name
  countdata_raw$dispersion = dlrs(x = countdata_raw$TN_ratio)
  countdata_raw
  
  #' segments obtained by Facets
  facets_segments = facets::procSample(x = countdata)
  facets_segments = facets_segments$out
  facets_segments = facets_segments[which(facets_segments$chrom == ROI[1, 'chromosome']), ]
  facets_segments
  
  output = list(countdata = countdata_raw, segments = facets_segments)
  output

}

#' Visualization of Facets output
plot_facets = function(data){
  data_raw = data$countdata
  data_cbs = data$segments
  
  #' ERBB2 coordinates:
  ERBB2_coordinates = data_raw$indx[which(data_raw$maploc >= 37842338 &
                                          data_raw$maploc <= 37886915)]
  #' raw log T/N distribution
  TN_raw = ggplot(data_raw, aes(x = indx, y = TN_ratio)) +
    geom_point() +
    theme(aspect.ratio = 0.5,
          axis.line.y = element_line(colour = 'black', size = 0.2),
          panel.background = element_blank()) +
    scale_y_continuous(limits = c(-2, 2)) +
    geom_hline(yintercept = seq(-2, 2, 1), size = 0.1, linetype = 'dashed') +
    labs(title = paste0('chromosome_17q', '; dlrs: ', round(data_raw$dispersion, 3)),
         y = 'Facets GC-corrected CnLR') +
    theme(plot.margin = margin(0, 0, 0, 0, "cm")) +
    geom_segment(aes(x = ERBB2_coordinates[1], xend = ERBB2_coordinates[length(ERBB2_coordinates)],
                     y = 1.3, yend = 1.3), color = 'red', size = 1.2) +
    geom_text(x = ERBB2_coordinates[1], y = 1.45, label = "ERBB2", vjust = 'middle')
  
  #' make the histogram
  TN_dispersion = ggplot(data_raw, aes(x = TN_ratio)) + 
    geom_histogram(aes(y = ..density..), bins = 50, colour="black", fill="white") +
    labs(title = paste0('median = ', round(median(data_raw$TN_ratio), 3), '; sd = ', 
                        round(sd(data_raw$TN_ratio), 3)),
         x = 'Facets-CnLR')
  
  
  #' CBS visualization:
  cbs_plot_raw = data_cbs
  name = substr(x = data_raw$sample, start = 17, stop = 33)
  
  genomdat = cbs_plot_raw$cnlr.median
  maploc = 1:length(genomdat)
  ii = cumsum(c(0, cbs_plot_raw$num.mark))
  mm = cbs_plot_raw$cnlr.median
  kk = length(ii)
  
  p2 = function(){
    par(
      mar = c(4, 2, 4, 2),
      mgp = c(2, 1, 0)
    )
    plot(1, 
         type = "n", 
         xlab = "",
         ylab = "", 
         xlim = c(0, max(ii)),
         ylim = c(-2, 2))
    segments(ii[-kk] + 1, 
             cbs_plot_raw$cnlr.median,
             x1 = ii[-1], 
             y1 = cbs_plot_raw$cnlr.median, 
             col = 'red')
    title(main = paste0('Facets-(default) segmentation; ', name), 
          xlab = 'index')
    abline(h = seq(-2, 2, 1), lty = 'dashed', lwd = 0.2)
  }
  
  #' make the output
  p1 = TN_raw + TN_dispersion
  plot_grid(p1, ggdraw(p2))
  
}

Facets_processed = readRDS('Data_out/Facets_segments.rds')

plot_list_facets = lapply(Facets_processed, function(x) plot_facets(x))



###############################################################################
#' Third panel: DryCleans Foreground (cleaned Tumors)
DryClean_metrics = function(data){
  file_name = basename(data)
  cleaned_data = as.data.frame(readRDS(data))
  cleaned_data$probes = paste(cleaned_data$seqnames, cleaned_data$start, sep = ';')
  cleaned_data = cleaned_data[which(cleaned_data$seqnames == ROI[1, 'chromosome']), ]
  cleaned_data = cleaned_data[which(cleaned_data$probes %in% RecognizedProbes$probes[which(RecognizedProbes$chromosome == ROI[1, 'chromosome'])]), ]
  cleaned_data = cleaned_data[which(cleaned_data$start > ROI[1, 'start'] &
                                        cleaned_data$start < ROI[1, 'end']), ]
  cleaned_data$TN_ratio = cleaned_data$foreground.log
  cleaned_data$indx = seq(from = 1, to = nrow(cleaned_data), by = 1)
  cleaned_data$sample = file_name
  cleaned_data$dispersion = dlrs(x = cleaned_data$TN_ratio)
  cleaned_data

}


rds = list.files('Tumor_cleaned/', full.names = T)
o = lapply(rds, function(x) DryClean_metrics(x))

#' DryClean segmentation
CBS_DryClean = function(data){
  data = as.data.frame(readRDS(data))
  name = basename(data)
  CBS_data = data[, c('seqnames', 'start', 'foreground.log')]
  CNA.object = CNA(genomdat = cbind(CBS_data$foreground.log),
                   chrom = CBS_data$seqnames, 
                   maploc = CBS_data$start, 
                   data.type = 'logratio', 
                   sampleid = name)
  smoothed_object = smooth.CNA(x = CNA.object)
  segmented_object = DNAcopy::segment(x = smoothed_object, 
                                      min.width = 5,
                                      undo.splits = 'sdundo', 
                                      undo.SD = 3)
  segmented_object
}


# lapply(rds, function(x) CBS_DryClean(x))




#' DryClean Visualization:
plot_dryclean = function(data_raw, data_cbs){
  file_name = basename(data_raw)
  data_raw = as.data.frame(readRDS(data_raw))
  data_raw$probes = paste(data_raw$seqnames, data_raw$start, sep = ';')
  data_raw = data_raw[which(data_raw$seqnames == ROI[1, 'chromosome']), ]
  data_raw = data_raw[which(data_raw$probes %in% RecognizedProbes$probes[which(RecognizedProbes$chromosome == ROI[1, 'chromosome'])]), ]
  data_raw = data_raw[which(data_raw$start > ROI[1, 'start'] &
                                        data_raw$start < ROI[1, 'end']), ]
  data_raw$indx = seq(from = 1, to = nrow(data_raw), by = 1)
  data_raw$sample = file_name
  data_raw$dispersion = dlrs(x = data_raw$foreground.log)

  #' raw log T/N distribution
  #' ERBB2 coordinates
  ERBB2_coordinates = data_raw$indx[which(data_raw$start >= 37842338 &
                                            data_raw$start <= 37886915)]
  TN_raw = ggplot(data_raw, aes(x = indx, y = foreground.log)) +
    geom_point() +
    theme(aspect.ratio = 0.5,
          axis.line.y = element_line(colour = 'black', size = 0.2),
          panel.background = element_blank()) +
    scale_y_continuous(limits = c(-2, 2)) +
    geom_hline(yintercept = seq(-2, 2, 1), size = 0.1, linetype = 'dashed') +
    labs(title = paste0('chromosome_17q', '; dlrs: ', round(data_raw$dispersion, 3))) +
    theme(plot.margin = margin(0, 0, 0, 0, "cm")) +
    geom_segment(aes(x = ERBB2_coordinates[1], xend = ERBB2_coordinates[length(ERBB2_coordinates)],
                     y = 1.3, yend = 1.3), color = 'red', size = 1.2) +
    geom_text(x = ERBB2_coordinates[1], y = 1.45, label = "ERBB2", vjust = 'middle')
  
  #' make the histogram
  TN_dispersion = ggplot(data_raw, aes(x = foreground.log)) + 
    geom_histogram(aes(y = ..density..), bins = 50, colour = "black", fill = "white") +
    labs(title = paste0('median = ', round(median(data_raw$foreground.log), 3), '; sd = ', 
                        round(sd(data_raw$foreground.log), 3)),
         x = 'foreground.log')
  
  
  #' concentrate on segmentation of cleaned tumor
  name = substr(x = file_name, start = 0, stop = 17)
  data_cbs = as.data.frame(readRDS(data_cbs))
  data_cbs$sample = file_name
  data_cbs = data_cbs[which(data_cbs$seqnames == ROI[1, 'chromosome']), ]
  
  genomdat = data_cbs$seg.mean
  maploc = 1:length(genomdat)
  ii = cumsum(c(0, data_cbs$num.mark))
  mm = data_cbs$seg.mean
  kk = length(ii)
  
  p2 = function(){
    par(
      mar = c(4, 2, 4, 2),
      mgp = c(2, 1, 0)
    )
    plot(1, 
         type = 'n',
         xlab = 'CBS_segments',
         ylab = 'seg.mean', 
         pch = '.', 
         ylim = c(-2, 2),
         xlim = c(0, max(ii)))
    
    abline(h = seq(-2, 2, 1), lty = 'dashed', lwd = 0.2)
    segments(ii[-kk] + 1,
             genomdat, 
             x1 = ii[-1], 
             y1 = genomdat, 
             col = 'red')
    
  title(main = paste0('CBS segmentation; ', name), 
          xlab = 'CBS_segments')
    abline(h = seq(-2, 2, 1), lty = 'dashed', lwd = 0.2)
  }
  
  #' make the output
  p1 = TN_raw + TN_dispersion
  plot_grid(p1, ggdraw(p2))
  
}

data_raw_paths = list.files('Tumor_cleaned/', all.files = F, full.names = T)
data_cbs_paths = list.files('Tumor_cleaned_CBS/', all.files = F, full.names = T)

plot_list_dryclean = list()
for(i in unique(data_raw_paths)){
  name = substr(i, start = 16, stop = 32)
  cbs = grep(pattern = name, x = data_cbs_paths, value = T)
  if(length(cbs) > 1){
    cbs = cbs[1]
  }
  plot_list_dryclean[[i]] = plot_dryclean(data_raw = i, data_cbs = cbs)
  
}


#' output the files
for(i in 1:65){
  ggsave_golden(plot_list_tn[[i]] / 
                  plot_list_facets[[i]] / 
                  plot_list_dryclean[[i]], filename = paste0('Figures/Sample_', i, '.pdf'), width = 12)
  
}

str(plot_list_tn[[1]])


#' out
