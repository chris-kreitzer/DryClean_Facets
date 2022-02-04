## Utility functions for DryClean processing
## 
## 10/08/2021
## chris-kreitzer
## 

dlrs = function(x) {                                                                                                                                                                                                                                        
  nx = length(x)                                                                                                                                                                                                                                          
  if (nx < 3) {                                                                                                                                                                                                                                              
    stop("Vector length>2 needed for computation")                                                                                                                                                                                                       
  }                                                                                                                                                                                                                                                        
  tmp = embed(x,2)                                                                                                                                                                                                                                        
  diffs = tmp[,2]-tmp[,1]                                                                                                                                                                                                                                 
  dlrs = IQR(diffs, na.rm = TRUE) / 1.34
  return(dlrs)                                                                                                                                                                                                                                             
}              


ggsave_golden = function(filename, plot, width, ...){
  ggsave(filename = filename, plot = plot, device = cairo_pdf, width = width, height = width / 1.61803398875)
}

CI_z <- function (x, ci = 0.95){
  `%>%` <- magrittr::`%>%`
  standard_deviation = sd(x)
  sample_size = length(x)
  Margin_Error = abs(qnorm((1-ci)/2))* standard_deviation/sqrt(sample_size)
  df_out = data.frame(sample_size = length(x), 
                      Mean = mean(x), 
                      sd=sd(x),
                      Margin_Error = Margin_Error,
                      'CI lower limit'= (mean(x) - Margin_Error),
                      'CI Upper limit'= (mean(x) + Margin_Error)) %>%
    tidyr::pivot_longer(names_to = "Measurements", 
                        values_to = "values", 1:6)
  return(df_out)
}
