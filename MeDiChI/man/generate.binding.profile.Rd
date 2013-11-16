\name{generate.binding.profile}
\alias{generate.binding.profile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Construct the binding profile ('kernel') required as input for
  'chip.deconv' and 'deconv.entire.genome'.}
\description{
  Construct a model-based binding profile as described in Reference (1).
}
\usage{
generate.binding.profile(fragment.distrib = function(x, ...) dgamma(x,
shape = 6, scale = 50), bs.size = 1, tile.size = 50, min.frag.size = 0,
positions = seq(0, 1001, by = 50), intensity.scaling = function(x, ...)
x, hybridization.prob = function(x, ...) as.integer(x > 10), interp = T,
plot = F, verbose = F, no.multicore=T, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fragment.distrib}{Distribution of DNA fragment lengths.}
  \item{bs.size}{Footprint of TF binding site on the genome.}
  \item{tile.size}{Length of the array probes (in bp).}
  \item{min.frag.size}{Cutoff for minimum DNA fragment size.}
  \item{positions}{The distances from the center of the profile for
    which to compute the relative intensity. A longer vector increases the
    running time of this procedure. Values for distances other than
    those listed here are computed via interpolation. See 'interp'.}
  \item{intensity.scaling}{Intensity of DNA fragments as a function of
    their length.}
  \item{hybridization.prob}{Hybridization probability as a function of
    length of complementary sequence.}
  \item{interp}{If TRUE, interpolate the values for distances other than
  those listed in 'positions'.}
  \item{plot}{If TRUE, plot the resulting profile.}
  \item{verbose}{If TRUE, be verbose.}
  \item{no.multicore}{Prevent use of multiple cores, even if 'multicore'
  is installed.}  
%  \item{kernel.method}{If 'mine' (the default), use the model described
%    in Reference (1); if 'jbd', use the method described in Reference (2). }
  \item{\dots}{Further parameters for 'fragment.distrib',
    'intensity.scaling', and 'hybridization.prob'.}
}
\details{
  No details.
}
\value{
  A two-column matrix containing positions (distance from center, in bp)
  and relative intensities of the profile.
}
\references{
  (1) Reiss, DJ and Facciotti, MT and Baliga, NS. (2007). "Model-based
  deconvolution of genome-wide DNA binding",
  Bioinformatics; doi: 10.1093/bioinformatics/btm592. \cr
  \url{http://baliga.systemsbiology.net/medichi}\cr
  (2) Qi, Y and et al. (2006). "High-resolution computational models of genome
  binding events", Nature Biotechnol, 24(8), 963-970. \url{http://cgs.csail.mit.edu/jbd}.
}
\author{
  David J Reiss, Institute for Systems Biology
  
  Maintainer: <dreiss@systemsbiology.org>
}
\seealso{
  chip.deconv, deconv.entire.genome, fit.peak.profile, generate.fake.data,
  MeDiChI-data
}
\examples{
## Compare profiles for DNA fragment distrs. with mean=300 and 400 bp.
kern.300 <- generate.binding.profile( fragment=function(x) dgamma( x,
                     shape=6, scale=50 ), verbose=TRUE )
kern.400 <- generate.binding.profile( fragment=function(x) dgamma( x,
                     shape=8, scale=50 ), verbose=TRUE )
plot( kern.300, typ="l", col="red" )
lines( kern.400, col="green" )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
