\name{MGSM27}
\alias{MGSM27}
\docType{data}
\title{
Myeloma Genome Signature Matrix 27
}
\description{
Newman et al's 2015 plus 5 myeloma specific cell types. Osteoclasts, Adipocytes, Osteoblasts, Multiple Myeloma Plasma Cells, and Plasma Memory Cells
}
\usage{data("MGSM27")}
\format{
  A data frame with 601 observations on the following 27 variables.
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
    \item{\code{MM.plasma.cell}}{a numeric vector}
    \item{\code{osteoblast}}{a numeric vector}
    \item{\code{osteoclast}}{a numeric vector}
    \item{\code{PlasmaMemory}}{a numeric vector}
    \item{\code{adipocyte}}{a numeric vector}
  }
}
\details{
MGSM27 as constructed for Identifying a High-risk Cellular Signature in the Multiple Myeloma Bone Marrow Microenvironment.
}
\source{
https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-3732/
https://www.ebi.ac.uk/arrayexpress/experiments/E-MEXP-3711/
https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-4152/
}
\examples{
data(MGSM27)
heatmap(as.matrix(MGSM27))
}
\keyword{datasets}
