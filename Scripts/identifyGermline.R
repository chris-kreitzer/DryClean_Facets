## Identify Germline events; Germline alterations
## 


#' normal.table.path must contain 
#' sample	
#' normal_cov	
#' decomposed_cov (from identify germline)
#' start_wash_cycle
grm = identify_germline(normal.table.path = "~/git/dryclean/inst/extdata/normal_table.rds", 
                        path.to.save = "~/git/dryclean/inst/extdata/", signal.thresh=0.5, pct.thresh=0.98)
