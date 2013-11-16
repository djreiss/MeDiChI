\name{generate.fake.data}
\alias{generate.fake.data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generate simulated ChIP-chip data.}
\description{
  Generate simulated ChIP-chip data, as described in Reference (1).
}
\usage{
generate.fake.data(posns = seq(1, 6001, by = 20), n.pts = 5, noise =
0.1, min.pk = 0.1, plot = F, verbose = F, in.pts = NULL, kernel = NULL,
tile.size = 100, reps = 1,
noise.func = function(data)noise*(1+sqrt(data)), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{posns}{Probe centers}
  \item{n.pts}{Number of (randomly-chosen) binding sites to add}
  \item{noise}{Level of noise (as described in Reference (1))}
  \item{min.pk}{Minimum peak intensity}
  \item{plot}{Plot the generated data?}
  \item{verbose}{Be verbose?}
  \item{in.pts}{2-column matrix with positions (column 1) and
    intensities (column 2) of input peaks; if n.pts is 'NA'.}
  \item{kernel}{Input peak profile to use.}
  \item{tile.size}{Input probe length (used to generate peak profile
    kernel if 'kernel' is 'NULL')}
  \item{noise.func}{Function used to generate noise as a function of
    signal intensity}
  \item{reps}{Number of replicate intensities per simulated probe}
  \item{\dots}{Further parameters for 'generate.binding.profile' (if
    'kernel' is 'NULL')}
}
\details{
  No details.
}
\value{
  A list of class containing the following three elements:
  \item{input}{A two-column matrix containing positions and intensities
    of input binding sites}
  \item{data}{A two-column matrix containing positions and intensities of simulated
    probes, can be passed directly to 'chip.deconv'}
  \item{kernel}{The peak profile kernel used to generate the data (as
    produced by 'generate.binding.profile')}
}
\references{
  Reiss, DJ and Facciotti, MT and Baliga, NS. (2007). "Model-based
  deconvolution of genome-wide DNA binding",
  Bioinformatics; doi: 10.1093/bioinformatics/btm592. \cr
  \url{http://baliga.systemsbiology.net/medichi}\cr
}
\author{
  David J Reiss, Institute for Systems Biology
  
  Maintainer: <dreiss@systemsbiology.org>
}
\seealso{chip.deconv, generate.binding.profile}
\examples{
## Generate data with 2 peaks at positions 3000 and 4000, deconvolve the
## data, and plot the resulting data and fit.
kern.300 <- generate.binding.profile( fragment=function(x) dgamma( x,
                     shape=6, scale=50 ), verbose=TRUE )
data <- generate.fake.data( in.pts=cbind( c( 3000, 4000 ), c( 1, 0.7 ) ), reps=3, 
                     kernel=kern.300 )
plot( data$data, pch=20 )
fit <- chip.deconv( data$data, center=NA, wind=NA, kernel=data$kernel,
                    fit.res=30, n.boot=10, verbose=TRUE, boot.sample.opt="case" )
plot( fit, boot="prob.scaled" )
print( fit )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
