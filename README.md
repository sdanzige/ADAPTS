# ADAPTS
Augments existing or de-novo cell-type signature matrices to deconvolve bulk gene expression data     This package expands on the techniques outlined in Newman et al.'s 2015 Nature Methods paper:      'Robust enumeration of cell subsets from tissue expression profiles'. to allow a user to easily add     their own cell types (e.g. a tumor specific cell type) to Newman's LM22 or other signature matrix.


To install this package in R, use devtools.

install.packages('devtools')

library(devtools)

devtools::install_github('sdanzige/ADAPTS')
