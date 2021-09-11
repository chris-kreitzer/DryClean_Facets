## Running the decompostion on tumor samples
## 





cov_out = start_wash_cycle(cov = coverage_file, 
                           detergent.pon.path = "~/git/dryclean/inst/extdata", 
                           whole_genome = TRUE, 
                           chr = NA, 
                           germline.filter = TRUE, 
                           germline.file = "~/git/dryclean/inst/extdata/germline.markers.rds")