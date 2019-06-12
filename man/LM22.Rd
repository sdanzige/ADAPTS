\name{LM22}
\alias{LM22}
\docType{data}
\title{
Leukocyte 22 data matrix
}
\description{
Newman et al.'s 2015 22 leukocyte signature matrix.
}
\usage{data("LM22")}
\format{
  A data frame with 547 observations on the following 22 variables.
  \describe{
    \item{\code{B.cells.naive}}{a numeric vector}
    \item{\code{B.cells.memory}}{a numeric vector}
    \item{\code{Plasma.cells}}{a numeric vector}
    \item{\code{T.cells.CD8}}{a numeric vector}
    \item{\code{T.cells.CD4.naive}}{a numeric vector}
    \item{\code{T.cells.CD4.memory.resting}}{a numeric vector}
    \item{\code{T.cells.CD4.memory.activated}}{a numeric vector}
    \item{\code{T.cells.follicular.helper}}{a numeric vector}
    \item{\code{T.cells.regulatory..Tregs.}}{a numeric vector}
    \item{\code{T.cells.gamma.delta}}{a numeric vector}
    \item{\code{NK.cells.resting}}{a numeric vector}
    \item{\code{NK.cells.activated}}{a numeric vector}
    \item{\code{Monocytes}}{a numeric vector}
    \item{\code{Macrophages.M0}}{a numeric vector}
    \item{\code{Macrophages.M1}}{a numeric vector}
    \item{\code{Macrophages.M2}}{a numeric vector}
    \item{\code{Dendritic.cells.resting}}{a numeric vector}
    \item{\code{Dendritic.cells.activated}}{a numeric vector}
    \item{\code{Mast.cells.resting}}{a numeric vector}
    \item{\code{Mast.cells.activated}}{a numeric vector}
    \item{\code{Eosinophils}}{a numeric vector}
    \item{\code{Neutrophils}}{a numeric vector}
  }
}
\source{
Newman, A. M. et al. Robust enumeration of cell subsets from tissue expression profiles. Nat. Methods 12, 453–457 (2015).
https://media.nature.com/original/nature-assets/nmeth/journal/v12/n5/extref/nmeth.3337-S2.xls
}
\examples{
data(LM22)
heatmap(as.matrix(LM22))
}
\keyword{datasets}
