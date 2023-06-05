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
XHMM: https://www.cell.com/action/showPdf?pii=S0002-9297%2812%2900417-X

## Early days of cytogenetic analysis:
Comparative Genomic Hybridization:    
- Comparative Genomic Hybridization for Molecular Cytogenetic Analysis of Solid Tumors (https://sci-hub.ru/https://www.science.org/doi/10.1126/science.1359641?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed)
https://mp.bmj.com/content/molpath/52/5/243.full.pdf  
- https://onlinelibrary.wiley.com/doi/epdf/10.1002/gcc.2870140405?saml_referrer    
- https://www.pnas.org/doi/epdf/10.1073/pnas.83.4.1031   
- CGH protocol: https://mp.bmj.com/content/molpath/52/5/243.full.pdf 

#### circular binary segmentation
- https://academic.oup.com/bioinformatics/article/27/15/2038/401729?login=false   
- https://watermark.silverchair.com/kxh008.pdf?token=AQECAHi208BE49Ooan9kkhW_Ercy7Dm3ZL_9Cf3qfKAc485ysgAAAtYwggLSBgkqhkiG9w0BBwagggLDMIICvwIBADCCArgGCSqGSIb3DQEHATAeBglghkgBZQMEAS4wEQQMvwwED3_wtBT8eN7MAgEQgIICiYGA9gTKbdM5JYjU46EpesVOQ_h4rUpy0LX6MdWGgESlOiykHL91QRkEV1D0UDngq8uQ39cn1rIOOxsBgtcTwBm-PwgOcsyU8wrS__W-WLrFm3g_iJMskUR7p4fpbXZuY3lkPt-gk8ebI1FHFSmr0krb8Zs_qxWlLUlOvsW5kTqcNNM_526v7aYmPoW_H_w4aZy87lAWbAhne0D1yLUj1ux3zhfju3g1zbEjlQS09-X5KMeoIuY17sxsNKyp2IvMhYid0q6tMLo-P1zkUMYPSqXQ43v4XvkhGr2SLeWYFq-g1E41CES1nSKN-2sI6M0qa7cDWYngcTedM1Lkjh9ybQ21982EzI5rczxPkPGo58T6PRiVEOKMXoHXYJzdOPSS1ScqHYwk1frin8_hku6d2ya6xQXsxAuw7er_wdX9dvu4fIqtGH71G157eUSZro1hJoEPmzuD5bieBxbI4LML42fw9-FHmKAdWZz8plnPCNHZccMDGP9kqoKdoxKbew7OHMbsvIAsmIK-vKQgOY2qtn5N0Qqa-3liim9mSGf19bHsLPJB5yxNX7N9M9Ds6Jf34cXUaAzC6Vcph9yuHPaOMR3kdpNSL6TCa7HOMJdGWZx-yaSlgcoJCnXRS4wXScPRF3PeicOFT9wY4jSSDP5PImO7ksgnLFNmb2oYgLlFobsQk-YAH_2ToskaGS0KsrBEEoU_43z_Fz94T8og5p_Wg4epGntggQgJuio-4wZOfEBYZJEZdfkB7vORlYoF5X22eFkLEkE1hBEqgb8mRWxOiN-N2Nd3-5vRP3ykPSyTJA8gVuKlAU3z1UWtqPmqcWj_PfaY58YoNCMA4craMj89J5j_0qBhKvZkYR0   
- https://www.cell.com/cms/10.1016/j.ajhg.2012.08.005/attachment/d70eab13-5bb1-44cf-ae79-83a33d359965/mmc1.pdf   



#### correcting CnLR for purity and ploidy: important ressource:
- https://www.nature.com/articles/s41586-020-2698-6.pdf   
- https://bitbucket.org/schwarzlab/refphase/src/master/R/calc_logr_thresholds.R
- https://github.com/lima1/PureCN/issues/40
- https://www.nature.com/articles/ng.2760#Sec8   

#### Purity-adjusted raw copy number:
- 2 * 2^(cnlr.median - dipLogR) - 2 * (1-purity)/purity    
- https://github.com/mskcc/facets/issues/7

#### FFPE sWGS (2014): important publication   
- https://pubmed.ncbi.nlm.nih.gov/25236618/   

#### TCGA copy number alterations - GISTIC analysis
- https://portals.broadinstitute.org/tcga/home


#### FACETS-preview:
https://bandla-chai.gitbook.io/facets-preview/installation2   


Marcin Imilienski   

Why do I want to be there? 
