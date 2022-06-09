# DryClean_Facets
A review on Facets and a potential implementation of DryClean

#### DryClean implementation:   
###### Creating PON
We start by creating a Panel of Normal (**PON**) for the DryClean pipeline;   
1,000 randomly sampled Breast Cancer; Breast Invasive Ductal Carcinoma; Gene Panel 468 (n = 2,705 samples) - retrieved from cBIO (07/07/2021) are used.    
Our sequencing approach (targeted panel; IMPACT-468) does not capture every genomic postion evenly throughout every sample; hence we need to perform several modifications.   
* We only keep positions where > 35 reads map to (in the normal sample; counts___merged as obtained from `snp-pileup`)   
* We only keep positions which are shared in at least 90% of **INPUT_NORMAL_SAMPLES**; this is necessary for creating the correct input_pon_matrix where n = genomic bins x m = samples; NOTE: the current algorithm considers positions (bins, probes) which are in at least 90% of input samples; if a sample does miss a specific position, this position will be supplemented with a '1'
* Probes coming from the Y-chromosome are not considered (only **1:22, X**)   
* We then perform a **mean normalization** to center the read counts (NOR.DP) around 0 (required for downstream analysis); Note: mean-normalization is based on a per-sample basis  
* The mean-normalization has the form of [x' = x / mean(x)]  

Modified dataframes are then converted into **GRanges objects**, where   
* start = actual position where sequencing read count was gathered
* end = start-position
* meta-column: **reads.corrected** == mean-normalized read counts   

Individual GRanges objects were stored locally, and subject to `prepare_detergent` with default parameters (`use.all = T, choose.randomly = F, choose.by.clustering = F`)


#### Update 03/31/2022:  
Marcin has created a messy genome (noise). Run dryclean on this messy simulation and see whether the background gets really flashed off. 

PCA to separate noise: Copy number variation detection and genotyping
from exome sequence data
Matrix decomposition: https://www.jstatsoft.org/article/view/v089i11
