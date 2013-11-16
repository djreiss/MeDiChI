\name{deconv.entire.genome}
\alias{deconv.entire.genome}
\alias{print.chip.deconv.entire.genome}
\alias{plot.chip.deconv.entire.genome}
\alias{coef.chip.deconv.entire.genome}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{High-resolution model-based deconvolution of normalized
  ChIP-chip data derived from tiling arrays.}
\description{
  Deconvolves an entire ChIP-chip data set, across all chromosomes,
  including multiple replicates, and running multiple bootstraps, if
  desired, by running 'chip.deconv' on multiple, contiguous, overlapping
  subsets of the data and combining the results.
}
\usage{
deconv.entire.genome(data, chroms=NA, window=6000, step.by=5500,
centers=NULL, quiet=F, plot=F, n.boot=1, fit.res=10, quant.cutoff="q0.85",
max.peak=NA, verbose=F, kernel=NA, no.multicore=T, ... )

print.chip.deconv.entire.genome(x, ...)

plot.chip.deconv.entire.genome(x, where = NA, center = NA, window = NULL, ...)

coef.chip.deconv.entire.genome(object, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Input data matrix, connection, or filename. See 'chip.deconv'.}
  \item{chroms}{Character vector listing chromosomes to include in
    analysis. The default ('NA') is to include all chromosomes. See
    'chip.deconv' for more information.}
  \item{window}{Size (in bp) of chunks to be processed
    one-at-a-time (see 'chip.deconv'). For
    'plot.chip.deconv.entire.genome', the coord. range that is plotted;
    the default ('NULL') is to use the 'window' that was input to
    'deconv.entire.genome()'. }
  \item{step.by}{Increment (in bp) in increment the 'center' coordinate
    of each 'chip.deconv' call.}
  \item{centers}{Optional vector of 'centers' to use (overrides
    'step.by' option. Vector may be optionally named by chromosome.}
%%  \item{fits.fin.only}{If TRUE, eturn only combined fits -- do not
%%    include the individual, 'raw' fits generated for each chunk.}
%%  \item{save.progress}{If non-NULL, periodically dump the fit to a file
%%    named by this parameter, enabling restart if the run fails at some
%%    point (rare, but this does happen!).}
%%  \item{in.fits}{If non-NULL, an input fit object, potentially saved in
%%    a prior run of this function (see 'save.progress') to re-start the
%%    deconvolution from.}
%%  \item{remove.near.edges}{If TRUE, remove any peaks near 5\% (in coord
%%    space) of chunk edges; their posns may not be correct.}
  \item{quiet}{If TRUE, be very quiet.}
  \item{plot}{If TRUE, plot each chunk as it is deconvolved.}
  \item{n.boot}{Number of bootstrap iterations. See 'chip.deconv'.}
  \item{fit.res}{Deconvolution resolution. See 'chip.deconv'.}
  \item{quant.cutoff}{See 'chip.deconv'.}
  \item{max.peak}{See 'chip.deconv'.}
  \item{verbose}{If TRUE, print out status messages. See 'chip.deconv'.}
  \item{kernel}{Required deconvolution kernel. See 'chip.deconv'.}
  \item{no.multicore}{Prevent use of multiple cores, even if 'multicore'
  is installed.}
  
  \item{x}{Object output from 'deconv.entire.genome()'}
  \item{where}{The chromosome(s) to be plotted. See Details.}
  \item{center}{Central chromosomal coordinate(s) to be plotted. See Details.}
  \item{object}{Object output from 'deconv.entire.genome()'}
  
  \item{\dots}{Further parameters passed to 'chip.deconv' or 'plot.chip.deconv'.}
}
\details{
  'deconv.entire.genome' breaks the full data set in to contiguous (in
  chromosomal coordinate) chunks of size 'window' and processes these
  using 'chip.deconv'. It increments the center of the 'window' by
  'step.by' base pairs, and combines the resulting fit by averaging
  nearby peaks in the overlapping segments. The defaults were chosen to
  contain overlaps of about 10\%, and constrain the required memory
  usage for each deconvolution.

  'deconv.entire.genome' distributes the processing of multiple sections of
  the chromosome to multiple processor cores, unless 'no.multicore' is
  set to TRUE.

  'plot.chip.deconv.entire.genome' will, by default, plot the entire fit
  across all chromosomes, via multiple plots, in chunks of size
  'window'. It is recommended to send this to a postscript or pdf device
  (rather than the screen). See 'plot.chip.deconv' for more information.
}
\value{
  A list of class 'chip.deconv.entire.genome', for which 'plot',
  'print', and 'coef'. This list is comprised of 'chip.deconv'
  objects (see 'chip.deconv') for each chromosome. See 'Example' for
  more information. The list is itemized as:
%%  \item{fits}{The individual raw fits to each chunk of each chromosome.
%%    Only included if 'fits.fin.only' is set to 'FALSE'.}
  \item{fits.fin}{The fits combined into a single one for each chromosome.}
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
\seealso{
  chip.deconv, plot.chip.deconv, fit.peak.profile, generate.fake.data,
  generate.binding.profile, MeDiChI-data,
  \code{lars}, \code{quadprog}, \code{Matrix}
}
\examples{
data( "halo.lowres", package="MeDiChI" )
## Deconvolve all three replicons: "Chr", "pNRC100", "pNRC200". Note
## this will take a while to run.

\dontrun{
fits <- deconv.entire.genome( data.halo.lowres, fit.res=30,
               n.boot=1, kernel=kernel.halo.lowres, verbose=FALSE, trace=FALSE )

## Plot the entire fit, across many windows:
plot( fits )

## Plot a piece of the finished product on just the chromosome:
plot( fits$fits.fin$Chr, center=20000, wind=10000 )
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
