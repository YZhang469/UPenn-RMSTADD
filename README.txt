Supplementary information / reproducible research files for the manuscript titled "Semiparametric Additive Modeling of the Restricted Mean Survival Time"

Authors: Y. Zhang and D. E. Schaubel
Code was written by Y. Zhang. Please address all questions or comments to yuan.zhang@pennmedicine.upenn.edu.

The code was written/evaluated in R with the following software information:
R version 4.2.2 (2022-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19044)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.utf8 
[2] LC_CTYPE=English_United States.utf8   
[3] LC_MONETARY=English_United States.utf8
[4] LC_NUMERIC=C                          
[5] LC_TIME=English_United States.utf8    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] survRM2_1.0-4  dplyr_1.0.10   ggplot2_3.4.0  haven_2.5.1    MASS_7.3-58.1 
[6] survival_3.4-0

loaded via a namespace (and not attached):
 [1] rstudioapi_0.14  knitr_1.40       magrittr_2.0.3   hms_1.1.2       
 [5] splines_4.2.2    tidyselect_1.2.0 munsell_0.5.0    colorspace_2.0-3
 [9] lattice_0.20-45  R6_2.5.1         rlang_1.0.6      fastmap_1.1.0   
[13] fansi_1.0.3      tools_4.2.2      grid_4.2.2       gtable_0.3.1    
[17] xfun_0.34        utf8_1.2.2       cli_3.4.1        withr_2.5.0     
[21] htmltools_0.5.3  ellipsis_0.3.2   yaml_2.3.6       digest_0.6.30   
[25] tibble_3.1.8     lifecycle_1.0.3  Matrix_1.5-1     vctrs_0.5.2     
[29] glue_1.6.2       evaluate_0.18    rmarkdown_2.18   compiler_4.2.2  
[33] pillar_1.8.1     generics_0.1.3   scales_1.2.1     forcats_0.5.2   
[37] pkgconfig_2.0.3 

Simulations (./simulation) were run in parallel on a Linux server with software information:
R version 4.2.2 (2022-10-31)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)

Matrix products: default
BLAS:   /misc/appl/R-4.2/lib64/R/lib/libRblas.so
LAPACK: /misc/appl/R-4.2/lib64/R/lib/libRlapack.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
 [9] LC_ADDRESS=C               LC_TELEPHONE=C
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base

other attached packages:
[1] survRM2_1.0-4  dplyr_1.0.10   ggplot2_3.4.0  haven_2.5.1    MASS_7.3-58.1
[6] survival_3.4-0

loaded via a namespace (and not attached):
 [1] magrittr_2.0.3   hms_1.1.2        splines_4.2.2    tidyselect_1.2.0
 [5] munsell_0.5.0    colorspace_2.1-0 lattice_0.20-45  R6_2.5.1
 [9] rlang_1.1.0      fansi_1.0.4      grid_4.2.2       gtable_0.3.1
[13] utf8_1.2.3       DBI_1.1.3        cli_3.5.0        withr_2.5.0
[17] ellipsis_0.3.2   assertthat_0.2.1 tibble_3.1.8     lifecycle_1.0.3
[21] Matrix_1.5-3     vctrs_0.5.1      glue_1.6.2       compiler_4.2.2
[25] pillar_1.8.1     forcats_0.5.2    generics_0.1.3   scales_1.2.1
[29] pkgconfig_2.0.3

This folder comprises three subfolders which contain the following files that can be used to reproduce all analysis results, tables, and figures of the manuscript:

./application/:
  A subfolder containing all scripts for data analysis reported in the manuscript (Section 5) and in the supplementary material (Appendix C).
  * biomj_postlt.sas7bdat: A simulated data set that has a similar structure to the liver transplant data analyzed in our manuscript.
  * analysis.r: An R script that performs the main data analysis, including generating the estimates of covariate effects with standard error and p-value (Table 3) and plotting the histogram of baseline mean RMST mu0's (Figures 1 and 2). To generate results in the supplementary material (Appendix C, Table C.1), change the parameter value L =  5*365 in the function "estBeta()" to L = 1*365 and L = 3*365 for 1-year RMST and 3-year RMST, respectively.
  * diagnostics.r: An R script that calculates C-statistic and Brier score of our proposed method and compare them with those of the additive RMST model, the stratified multiplicative RMST model and the Cox PH model (Table 4).

./function/:
  A subfolder containing a file "estBeta.r" which calls the estimating function "estBeta()". The files in both "./application" and "./simulation" folders source this file.

./results/:
  A subfolder containing all simulation and data analysis results.

./simulation/:
  A subfolder containing all scripts for simulation studies reported in the manuscript (Section 4) and in the supplementary material (Appendix B).
  * sim.r: An R script that performs the main simulation study (Table 1).
  * sim_sensitivity.r: An R script that performs the sensitivity analysis (Table 2).
  * sim_supp.r: An R script that performs the simulation study with a considerably large number of covariates (Appendix B, Table B.1).
