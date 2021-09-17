## Identify germline alteration from the normals;
## 
## working with default parameters
## input is a table which contains the normals as well as the 
## decomposed normals from start_wash_cycle()
## importantly - normals are treated like tumors

## 09/17/2021
## chris-kreitzer


library(parallel)
library(dryclean)


## the input data table must be in .rds format
## moreover it must be a data.table::data.table() object
## columns that are required
#' sample
#' normal_cov
#' decomposed_cov

working_path = ' '

Germline = identify_germline(normal.table.path = paste0(working_path, 'decomposed_samples/normal_table.rds'), 
                             path.to.save = paste0(working_path, 'decomposed_samples'), 
                             signal.thresh = 0.5, 
                             pct.thresh = 0.98)
