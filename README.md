# DryClean_Facets
A review on Facets and a potential implementation of DryClean

#### DryClean implementation:   
We start by creating a Panel of Normal (**PON**) for the DryClean pipeline;   
99 randomly sampled Breast Cancer; Breast Invasive Ductal Carcinoma; Gene Panel 468 (n = 2,705 samples) - retrieved from cBIO (07/07/2021) are used.    
Our sequencing approach (targeted panel; IMPACT-468) does not capture every genomic postion evenly throughout every sample; hence we need to perform several modifications.   
* We only keep positions where > 35 reads map to (in the normal sample; counts___merged as obtained from `snp-pileup`)   
* We only keep positions which are shared among all **INPUT_NORMAL_SAMPLES**; this is necessary for creating the correct input_pon_matrix where n = genomic bins x m = samples   
* We then perform a **mean normalization** to center the read counts (NOR.DP) around 0 (required for downstream analysis)  
* The mean-normalization has the form of [x' = x - mu / (max(x) - min(x)]   
* 
*  
 
