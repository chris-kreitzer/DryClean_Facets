## The second step in DryClean is to identify germline events,
## in order to remove them from the panel of normals;



## 09/06/2021
## chris-kreitzer





decomp.1 = start_wash_cycle(cov = sample.1, detergent.pon.path = "~/git/dryclean/inst/extdata/", whole_genome = TRUE, chr = NA, germline.filter = FALSE)