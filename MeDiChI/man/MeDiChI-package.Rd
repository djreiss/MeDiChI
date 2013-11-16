\name{MeDiChI-package}
\alias{MeDiChI-package}
\alias{MeDiChI}
\docType{package}
\title{
MeDiChI ChIP-chip deconvolution library
}
\description{
Deconvolution of Chromatin-IP-to-microarray (ChIP-chip) data to identify
binding sites at high resolution across the genome.\cr
To be published in Bioinformatics and submitted to Bioconductor.\cr
This package includes all data described in the manuscript.
}
\details{
\tabular{ll}{
Package: \tab MeDiChI\cr
Type: \tab Package\cr
Version: \tab 0.1-4\cr
Date: \tab 2007-11-25\cr
License: \tab GPL version 3\cr
}
demo(MeDiChI)\cr
}
\author{
  David J Reiss, Institute for Systems Biology
  
  Maintainer: <dreiss@systemsbiology.org>
}
\references{
  (1). Reiss, DJ and Facciotti, MT and Baliga, NS. (2007). "Model-based
  deconvolution of genome-wide DNA binding",
  Bioinformatics; doi: 10.1093/bioinformatics/btm592. \url{http://baliga.systemsbiology.net/medichi}\cr
  (2) Qi, Y and et al. (2006). "High-resolution computational models of genome
  binding events", Nature Biotechnol, 24(8), 963-970. \url{http://cgs.csail.mit.edu/jbd}.
}
\keyword{ package }
\seealso{
  chip.deconv, deconv.entire.genome, fit.peak.profile, generate.fake.data,
  generate.binding.profile, MeDiChI-data,
  \code{<lars>}, \code{<quadprog>}, \code{<Matrix>}
}
\examples{
demo( MeDiChI )

## Run the demo yourself:
data( "halo.lowres", package="MeDiChI" )

fit <- chip.deconv( data.halo.lowres, where="Chr", fit.res=30,
               center=650000, wind=20000, max.steps=100, n.boot=10,
               kernel=kernel.halo.lowres, verbose=TRUE, boot.sample.opt="case" )

coef( fit ) ## Print out the coefficients
plot( fit, plot.genes=TRUE, cex=0.5, cex.lab=0.8, cex.axis=0.8 ) 
}
