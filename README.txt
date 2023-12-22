Supplementary information / reproducible research files for
Title: "A Bayesian Basket Trial Design Using Local Power Prior"

Authors: Haiming Zhou, Rex Shen, Sutan Wu and Philip He
Code was written by Haiming Zhou and Rex Shen. 
In case of questions or comments please contact haiming2019@gmail.com

The code was written/evaluated in R with the following software versions:
R version 4.2.2 (2022-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19045)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8   
[3] LC_MONETARY=English_United States.utf8 LC_NUMERIC=C                          
[5] LC_TIME=English_United States.utf8    

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] cluster_2.1.4     plyr_1.8.7        partitions_1.10-7 basket_0.10.11    bhmbasket_0.9.5  
 [6] R2jags_0.7-1      rjags_4-13        coda_0.19-4       doParallel_1.0.17 iterators_1.0.14 
[11] foreach_1.5.2    

loaded via a namespace (and not attached):
 [1] tidyselect_1.2.0   purrr_0.3.4        graphlayouts_1.0.0 lattice_0.20-45    colorspace_2.0-3  
 [6] vctrs_0.5.2        generics_0.1.3     viridisLite_0.4.0  gmp_0.6-5          utf8_1.2.2        
[11] rlang_1.1.1        pillar_1.8.0       glue_1.6.2         withr_2.5.0        tweenr_2.0.2      
[16] RColorBrewer_1.1-3 lifecycle_1.0.3    munsell_0.5.0      gtable_0.3.0       codetools_0.2-18  
[21] fansi_1.0.3        itertools_0.1-3    tidygraph_1.2.3    Rcpp_1.0.10        polynom_1.4-1     
[26] scales_1.2.0       abind_1.4-5        farver_2.1.1       gridExtra_2.3      ggforce_0.4.1     
[31] ggplot2_3.4.3      digest_0.6.29      dplyr_1.1.0        ggrepel_0.9.1      rbibutils_2.2.8   
[36] mathjaxr_1.6-0     polyclip_1.10-4    grid_4.2.2         Rdpack_2.4         cli_3.6.0         
[41] tools_4.2.2        magrittr_2.0.3     tibble_3.1.7       ggraph_2.1.0       crayon_1.5.1      
[46] tidyr_1.2.0        pkgconfig_2.0.3    GenSA_1.1.8        ellipsis_0.3.2     MASS_7.3-58.1     
[51] rstudioapi_0.13    viridis_0.6.2      R6_2.5.1           boot_1.3-28        R2WinBUGS_2.1-21  
[56] igraph_1.3.4       compiler_4.2.2    


This folder contains the following data and files that can be used to reproduce all analysis results of the manuscript.
It contains four subfolders containing the following files:

./case_study/:
    code_case_study_alpha=0.05.R
    An R script that contains the code of the analysis reported in the paper (section 4) when type I error is controlled at 0.05. The resulted workspace is saved in case_study_alpha=0.05.RData.

    code_case_study_alpha=0.1.R
    An R script that contains the code of the analysis reported in the paper (section 4) when type I error is controlled at 0.1. The resulted workspace is saved in case_study_alpha=0.1.RData.  

./functions/:
    functions.R
    An R script that contains all functions to implement each Bayesian information borrowing method on a single dataset.

    cluster.R
    An R script that contains functions required for local MEM method, copied from Liu et al (2022).

    utils.R
    An R script that contains functions required for BCHM, copied from Chen and Lee (2020).

    functions_sim.R
    An R script that contains functions for simulations.

    functions_parallel.R
    An R script that contains functions for parallel computing in simulations.

./misc/:
    localPP_SimilarityMatrix.R
    An R script that was used to generate the similarity matrix of local PP in the manuscript (section 2.1).
    
./simulation/
    code_sim.R
    An R script that performs the simulations reported in the manuscript (section 3). Simulation was performed in parallel on 5 cores on Windows 10. Results of the simulation are saved into the folder ./intermediate_results/ since computation of the whole simulation can take 2 days since some MCMC methods are time-consuming.
    For faster computation the number of replications within each simulation scenario could be reduced from nperclust=1000. If the machine used for computation is not able to compute on 5 cores simultaneously the number of cores used for simulation should be reduced by changing the nclust argument.
    
    ./intermediate_results/
    A folder containing the results of code_sim.R
    
    results_sim.R
    An R script that takes the RData files from the intermediate_results and creates the
    results tables presented in the manuscript.
    
