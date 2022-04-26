## Investiage the detergent (L and S matrix)
## 
## start: 04/26/2022
## chris-kreitzer

clean()
gc()
library(RColorBrewer)
library(colorspace)
library(circlize)

detergent = readRDS('detergent.rds')
PON = vroom::vroom('~/Documents/MSKCC/07_FacetsReview/DryClean/DataProcessed/PON_normalized.txt', delim = '\t')

L = detergent$L
S = detergent$S
template = detergent$template
template = as.data.frame(template)
template$bin = paste(template$seqnames, template$start, sep = ';')
germline = detergent$inf_germ
germline = as.data.frame(germline)
germline_events = germline[which(germline$germline.status == T), ]
germline_events$bin = paste(germline_events$seqnames, germline_events$start, sep = ';')
germline_events = germline_events[which(germline_events$black_list_pct > 0.79), ]
dim(germline_events)

row.names(L) = template$bin
row.names(S) = template$bin


## how does Germline events look like:
germline_l = L[which(row.names(L) %in% germline_events$bin), ]
germline_s = S[which(row.names(S) %in% germline_events$bin), ]
germline_p = PON[which(PON$duplication %in% c('1;11168048', '1;39321150', '22;41545611')), ]

colnames(germline_l) = colnames(germline_p)[2:ncol(germline_p)]

colr = colorRamp2(
  seq(min(germline_l), max(germline_l), length = 3),
  c("#000099", "#EEEEEE", "#FF0000"),
  space = "RGB")

ComplexHeatmap::Heatmap(germline_s,
                        col = colr,
                        cluster_rows = T, 
                        cluster_columns = F, 
                        row_names_gp = gpar(fontsize = 6), 
                        show_column_names = F,
                        show_row_names = F, 
                        name = 'germline')



dim(germline_s)
xx = PON[which(PON$duplication %in% c('1;17084300', '2;189541550', '13;42829970')), ]
View(xx)
detergent$inf_germ

