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
