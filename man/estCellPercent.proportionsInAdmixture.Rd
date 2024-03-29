% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/onlyDeconAlgorithms.R
\name{estCellPercent.proportionsInAdmixture}
\alias{estCellPercent.proportionsInAdmixture}
\title{WGCNA::proportionsInAdmixture deconvolution}
\usage{
estCellPercent.proportionsInAdmixture(refExpr, geneExpr, marker_set = NULL)
}
\arguments{
\item{refExpr}{a data frame representing immune cell expression profiles. Each row represents an expression of a gene, and each column represents a different immune cell type. colnames contains the name of each immune cell type and the rownames includes the genes' symbol. The names of each immune cell type and the symbol of each gene should be unique. Any gene with missing expression values must be excluded.}

\item{geneExpr}{a data frame representing RNA-seq or microarray gene-expression profiles of a given complex tissue. Each row represents an expression of a gene, and each column represents a different experimental sample. colnames contain the name of each sample and rownames includes the genes' symbol. The name of each individual sample and the symbol of each gene should be unique. Any gene with missing expression values should be excluded.}

\item{marker_set}{data frames of one column, that includes a preselected list of genes that likely discriminate well between the immune-cell types given in the reference data. (DEFAULT: NULL, i.e. one for each gene in the refExpr)}
}
\value{
A matrix with cell type estimates for each samples
}
\description{
Use R function proportionsInAdmixture to estimate the cell count percentage
Uses the 'WGCNA' package

cellEst <- estCellPercent.proportionsInAdmixture(refExpr)
}
\examples{
\donttest{
#This toy example 
library(ADAPTS)
fullLM22 <- ADAPTS::LM22[1:30, 1:4]
smallLM22 <- fullLM22[1:25,] 

cellEst <- estCellPercent.proportionsInAdmixture(refExpr=smallLM22, geneExpr=fullLM22)
}
}
