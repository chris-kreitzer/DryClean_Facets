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