## Identify germline alteration from the normals;
## 
## 


#' normal.table.path must contain 
#' sample	
#' normal_cov	
#' decomposed_cov (from identify germline)
#' start_wash_cycle




Germline = identify_germline(normal.table.path = paste0(working_path, 'decomposed_samples/normal_table.rds'), 
                             path.to.save = paste0(working_path, 'decomposed_samples'), 
                             signal.thresh = 0.5, 
                             pct.thresh = 0.98)
