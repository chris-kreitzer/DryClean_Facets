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


#' TN metrics;
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


files = list.files('Tumor_countsFile/', full.names = T, all.files = F, include.dirs = T)
files = paste0('/Users/chriskreitzer/Documents/GitHub/DryClean_Facets/', files)


x = lapply(files, function(x) TN_metrics(x))
y = data.table::rbindlist(x)


#' Segmentation; CBS
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




#' Visualizations and noise quantifications
plot_samples = function(data_raw, data_cbs){
  TN_raw = ggplot(data_raw, aes(x = indx, y = TN_ratio)) +
    geom_point() +
    theme(aspect.ratio = 0.5,
          axis.line.y = element_line(colour = 'black', size = 0.2),
          panel.background = element_blank()) +
    scale_y_continuous(limits = c(-2, 2)) +
    geom_hline(yintercept = seq(-2, 2, 1), size = 0.1, linetype = 'dashed') +
    labs(title = 'chromosome_17q') +
    theme(plot.margin = margin(0, 0, 0, 0, "cm"))
  TN_dispersion = ggplot(data_raw, aes(x = TN_ratio)) + 
    geom_histogram(aes(y = ..density..), bins = 50, colour="black", fill="white")
  TN_raw + TN_dispersion
    # geom_segment(aes(x = ERBB2_coordinates[1], xend = ERBB2_coordinates[length(ERBB2_coordinates)],
    #                  y = 1.1, yend = 1.1), color = 'red', size = 1.2) +
    # geom_text(x = ERBB2_coordinates[1], y = 1.15, label = "ERBB2", vjust = 'middle')
  
}

DNAcopy::
plot_list = list()
for(i in unique(y$sample)){
  data_raw = y[which(y$sample == i), ]
  plot_list[[i]] = plot_samples(data_raw = data_raw)
}

library(cowplot)

plot_grid(plot_list[[1]], ggdraw(p1))



p1 = function() {
  par(
    mar = c(8, 2, 8, 2),
    mgp = c(2, 1, 0)
  )
  plot(a[[1]])
}

ggdraw(p1)








plot(a[[1]])





plot_list$`countsMerged____P-0000584-T03-IM6_P-0000584-N01-IM6.dat.gz`
plot_list$`countsMerged____P-0003195-T02-IM6_P-0003195-N01-IM6.dat.gz`

a = lapply(x, function(x) CBS_segmentation(x))

plot(a[[2]])
plot(x = a[[1]], 
     xmaploc = F, altcol = TRUE, sbyc.layout= NULL, 
     pt.pch = '.', pt.cex  = 2, 
     pt.cols = c('blue', 'grey'),
     segcol = 'brown', ylim = c(-2, 2),
     xlab = 'index chr18', main = '')

abline(h = seq(-2, 2, 1), lty = 'dashed', lwd = 0.5)
title(xlab = 'INDEX chr17')
title(main = "Main title", sub = "Sub-title",
      xlab = "X axis", ylab = "Y axis",
      cex.main = 2,   font.main= 4, col.main= "red",
      cex.sub = 0.75, font.sub = 3, col.sub = "green",
      col.lab ="darkblue"
)

o = readRDS('~/Desktop/re.rds')
plot(o[[1]])


lapply(x, function(x) CBS_segmentation(data = x))

sample_CBS = list()
for(i in unique(y$sample)){
  sample_CBS[[i]] = CBS_segmentation(data = y[which(y$sample == i), ])
}



plot(sample_CBS[[2]])
