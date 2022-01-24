# Improving Adaptive Seamless Designs through Bayesian optimization

Paper on arXiv: https://arxiv.org/abs/2105.09223

# README

This repo contains the code to generate the plots, tables for this paper, the paper pdf itself as well as the code to run the simulations.
Since running the simulations needs quite a few resources the results of the simulations are also contained.

## Files

- `./arxiv`: Arxiv version of the paper.
- `./benchmark`: Contains the code that runs the benchmarks (`batchtools_grid.R`) and creates the results given in `./benchmark_results`. `setwd('./benchmark')` to run the code.
- `./benchmark_results`: Results of the benchmark experiments.
- `./generator`: Code that creates the tables and figures from the files in `./benchmark_results`. `setwd('.')` (the directory of this README) to run the code.
- `./generated`: Figures and tables created by the code in `./generator`.
- `Makefile`:
	- `make plots` generates plots and tables.
	- `make pdf` generates the `main.pdf`.
	- `make all` runs the two above.
	- Use `make tinytex` in case some tex packages are missing. Tinytex can deal with that and generate the `main.pdf` as well.
	- `make benchmark` runs the expensive experiments. Only recommended if you have set up [`batchtools`](https://mllg.github.io/batchtools/) with access to a cluster.

## Packages

To run the _generators_ you have to install:

```r
pkgs = c("data.table","ggplot2","scales","stringi","knitr","kableExtra","dplyr","mlr3misc","gridExtra","memoise","smoof","batchtools","FNN","tidyr","colorspace","scales","BBmisc")
install.packages(pkgs)
```

Also pdfcrop is needed:
```
sudo apt-get install texlive-extra-utils  
```
```
pacman -S texlive-core
```

## Session Info

```
R version 3.6.3 (2020-02-29)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 21.04

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=de_DE.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=de_DE.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] grid      stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] smoof_1.6.0.2     checkmate_2.0.0   ParamHelpers_1.14 memoise_2.0.1    
 [5] gridExtra_2.3     mlr3misc_0.10.0   dplyr_1.0.7       kableExtra_1.3.4 
 [9] knitr_1.37        stringi_1.7.6     scales_1.1.1      ggplot2_3.3.5    
[13] data.table_1.14.2

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.7         svglite_2.0.0      FNN_1.1.3          tidyr_1.1.4       
 [5] prettyunits_1.1.1  digest_0.6.28      utf8_1.2.2         R6_2.5.1          
 [9] backports_1.3.0    evaluate_0.14      httr_1.4.2         pillar_1.6.4      
[13] progress_1.2.2     rlang_0.4.12       lazyeval_0.2.2     misc3d_0.9-1      
[17] rstudioapi_0.13    rmarkdown_2.10     labeling_0.4.2     mco_1.15.6        
[21] webshot_0.5.2      stringr_1.4.0      htmlwidgets_1.5.3  munsell_0.5.0     
[25] compiler_3.6.3     xfun_0.29          pkgconfig_2.0.3    systemfonts_1.0.2 
[29] BBmisc_1.11        htmltools_0.5.1.1  tcltk_3.6.3        tidyselect_1.1.1  
[33] tibble_3.1.6       batchtools_0.9.15  fansi_1.0.2        viridisLite_0.4.0 
[37] crayon_1.4.2       withr_2.4.3        rappdirs_0.3.3     jsonlite_1.7.3    
[41] gtable_0.3.0       lifecycle_1.0.1    magrittr_2.0.1     cachem_1.0.5      
[45] farver_2.1.0       xml2_1.3.2         brew_1.0-6         ellipsis_0.3.2    
[49] generics_0.1.0     vctrs_0.3.8        fastmatch_1.1-3    plot3D_1.4        
[53] RColorBrewer_1.1-2 tools_3.6.3        RJSONIO_1.3-1.4    glue_1.6.0        
[57] purrr_0.3.4        hms_1.1.1          parallel_3.6.3     fastmap_1.1.0     
[61] colorspace_2.0-2   base64url_1.4      rvest_1.0.1        plotly_4.10.0     
```