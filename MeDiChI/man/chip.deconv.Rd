\name{chip.deconv}
\alias{chip.deconv}
\alias{print.chip.deconv}
\alias{plot.chip.deconv}
\alias{coef.chip.deconv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{High-resolution model-based deconvolution of normalized
  ChIP-chip data derived from tiling arrays.}
\description{
  Deconvolves a subset of data on one chromosome, including multiple
  replicates, and running multiple bootstraps, if desired. To deconvolve
  an entire data set across multiple chromosomes, see 'deconv.entire.genome'.
}
\usage{
chip.deconv(data, where = NA, center = NA, 
window = 30000, fit.res = 10, max.steps = 200, 
post.proc.factor = 2, min.npeaks = 0, max.npeaks = 99999,
selection.method = "bic", quant.cutoff = "q0.85", 
n.boot = 1, boot.sample.opt = c("residual","resample","case","wild",
"position","replicate")[1], max.peak = NA, boot.vary.res = F,
kernel = NA, tile.distance = NA, verbose = T, trace = F, ...)
%% cluster.nodes = NA, plot.status = F, ls.final.do = T, interp = T, plot.boot = F, 

plot.chip.deconv(x, boot.results = c("scaled.prob", "prob", "scale",
"conf=95", "NONE")[1], where = NA, center = NA, window = NULL, verbose = F,
plot.genes = F, org = NA, hi.res = NA, quants = c(0.95, 0.5, 0.05), 
smooth = T, ylim=NULL, ... )
%%pch=20, cex=0.5, ...)
%%scale.boot.results = NA, 

print.chip.deconv(x, ...)

coef.chip.deconv(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Input data matrix, connection, or filename. See Details.}
  \item{center}{Central chromosomal coordinate for the subset of data to
    deconvolve, or to be plotted. See Details.}
  \item{where}{The chromosome of the subset of data to deconvolve. See
    Details.}
  \item{window}{The window size (in base-pairs) of the subset of
    data to deconvolve; see Details. For 'plot.chip.deconv', limit the
    window size that is plotted. The default ('NULL'), is to use the
    'window that was input to 'chip.deconv()'. }
  \item{fit.res}{Desired deconvolution resolution (base-pairs).}
  \item{kernel}{Required deconvolution kernel. See Details.}
  \item{max.steps}{Limit the number of LARS steps taken. Should opt
    towards higher values (above 200), as the fit from the highest step
    is used to estimate the noise for BIC.}
  \item{post.proc.factor}{Post-processing filter for combining
    deconvolution coefficients, in units of (n*fit.res).}
  \item{min.npeaks}{Optionally limit the minimum number of
    coefficients. Default is '0' -- no lower limit.}
  \item{max.npeaks}{Optionally limit the maximum number of
    coefficients. Default is '99999' -- no upper limit.}
  \item{selection.method}{Use argmin(method) to choose optimal model, 
    one of c('bic','aic'). Default is 'bic'.}
  \item{quant.cutoff}{Intensity or quantile cutoff for data to be processed. Limits
    the locations of potential sites to those near probes that are above
    this quantile. See Details. Default is 'q0.85' -- use 85th quantile.}
  %%\item{ls.final.do}{If TRUE, perform final positive-only least-squares
  %%  optimization on parameters for more accurate estimates of peak
  %%  intensities.}
  \item{n.boot}{Number of bootstrap iterations to perform. If '0'
    (default) then no bootstraps are performed.}
  \item{boot.sample.opt}{Bootstrap resampling option. See
    Details. Default is 'case.}
  \item{max.peak}{Any coefficient with an intensity above this
    threshold is set to this value. Default is 'NA' -- no cutoff. }
  \item{boot.vary.res}{If TRUE, vary the resolution (around 'fit.res') during
    bootstraps. May result in more realistic solutions.}
%%  \item{cluster.nodes}{Optionally distribute boostrap runs across
%%    processors or cluster nodes. See Details.}
%%  \item{interp}{If TRUE, interpolate probe intensities for data where
%%    'tile.distance' is >= 300bp (this is necessary to prevent unrealistic
%%    solutions).}
  \item{tile.distance}{Distance (in base pairs) between adjacent probes
    on the array. If 'NA' (the default) then this is computed from the data.}
%%  \item{plot.boot}{If TRUE, plot each bootstrap fit as it is
%%    performed.}
  %%\item{plot.status}{If TRUE, generate plots as the fit is performed,
  %%  e.g. including data that is included in fit (above
  %%  'quantile.cutoff', and BIC and AIC curves).}
  \item{verbose}{If TRUE, print out status messages.}
  \item{trace}{If TRUE, print out LARS progress.}
%%  \item{no.multicore}{If TRUE, prevent parallelizatin of bootstraps to
%%    multiple processor cores.}

  \item{x}{Object output from 'chip.deconv()'}
  \item{object}{Object output from 'chip.deconv()'}
  \item{boot.results}{Plot bootstrap distributions, either 'prob'
    (posterior probability), 'scaled.prob' (intensity-scaled posterior
    probability), 'conf=95' (95\% confidence intervals), or 'NONE' (no
    bootstrap results plotted). Default is 'scaled.prob'.}
%%  \item{center}{Limit the central coordinate of the window that is
%%    plotted. If 'NA' (the default), then the center that was input to
%%    'chip.deconv()' is used.}
%%  \item{verbose}{If TRUE, be verbose.}
  \item{plot.genes}{If TRUE, include gene positions in the
    plot. Currently only supported for Halobacterium and S. cerevisiae.}
  \item{org}{Organism used for 'plot.genes', one of either 'halo' or 'yeast'.}
  \item{quants}{Quantiles of model fits from bootstraps to include in
    plot (see Figure 3 in paper). Default is to plot 5th, 50th, and 95th
    quantiles.}
  \item{hi.res}{Compute high-resolution model fits at this
    resolution. Default is 'fit.res' as provided to 'chip.deconv()'.}
  \item{smooth}{If TRUE, plot kernel-smoothed bootstrap distributions, rather
    than simple counts.}
  \item{ylim}{Adjust the y-limits of the plot.}
%%  \item{pch}{Plot character. See 'par'.}
%%  \item{cex}{Plot character size. See 'par'.}
  \item{\dots}{Additional parameters passed to 'lars' or 'plot' or
    'read.table' (if the 'data' parameter is a file name).}
}
\details{
  'chip.deconv' is used to deconvolve a limited subset of the entire
  tiling array data set. Usually this subset is limited to a single
  chromosome (the 'where' parameter) and a range of coordinates on the
  chromosome between 'center - window / 2' and 'center + window /
  2'. The method will identify potential binding sites (coefficients)
  across this range at a resolution given by 'fit.res', up to a maximum
  of n.coeffs = window / fit.res coefficients.
  
  It is not recommended (due to memory constraints) to attempt to fit
  more than a few thousand potential coefficients using this
  function. This number may be decreased by either decreasing the
  'window' size or increasing the 'fit.res'. To deconvolve larger
  windows at high resolution, use 'deconv.entire.genome'.
  
  The input 'data' can be formatted as one of:\cr
  (a) a 2-column matrix or data-frame containing probe coordinates
  (column 1) and intensities (column 2), with optional rownames containing
  the chromosome identifier for each probe, or\cr
  (b) a 3-column data frame containing probe chromosome identifiers
  (column 1), coordinates (column 2) and intensities (column 3), or\cr
  (c) a file name or connection pointing to either a GFF file or a
  3-column tab-delimited file formatted as in (b).\cr

  Replicates are included in the input data simply as additional rows
  with the same probe coordinate.\cr

  Note that 'intensities' in the above description refers to relative
  (potentially normalized) intensities, exponentiated log-ratios
  (e.g. versus a reference) or some such quantity. The data should not
  be logged, as this will not conform with the MeDiChI peak model (see
  the reference below for details).\cr
  
  If the entire input 'data' matrix is to be processed using this
  function (only recommended if a subset of the entire array has been
  pre-selected; see above), or if the data only cover one chromosome,
  then the chromosome identifiers are not required.\cr
  
  If the input 'data' contains the entire array (including multiple
  chromosomes) then only the contiguous subset of probes within the
  window described above and with the matching chromosome identifier
  will be deconvolved.\cr

  'quant.cutoff' can be either a character starting with "q",
  e.g. "q0.85" to represent a quantile cutoff, or a numeric value,
  representing an absolute intensity cutoff. Only positions within +/-
  one 'tile.size' of any probe that has an intensity greater than this
  cutoff are considered to contain potential binding sites. This
  parameter may be used to decrease the runtime of the function.\cr

  'boot.sample.opt' can be one of 'wild' (Default) -- wild resampling (see
  \url{http://en.wikipedia.org/wiki/Bootstrapping_(statistics)\#Wild_bootstrap});
  'residual' -- wild resampling, only run on residuals for estimation of
  coefficient p-values; 'case' -- case resampling; 'position' -- resample
  the central positions of the probes; 'replicate' -- resample probe
  intensities from the range given by their replicates (if replicates
  exist); or 'resample' -- case sampling of intensities only for
  estimation of coefficient p-values.\cr

  'kernel' is a 2-column matrix providing position (column 1) and
  intensity (column 2) of the deconvolution kernel (profile model) to be
  used. This may be computed using 'generate.binding.profile' and
  parameters for this model may be learned from the data using
  'fit.peak.profile'.\cr

%%  'chip.deconv distributes the processing of multiple bootstraps
%%  to multiple processor cores, unless 'no.multicore' is set to FALSE.
  
%%  'cluster.nodes' is a list of cluster node names to distribute the
%%  bootstrap runs on. For example to distribute on a single quad-core
%%  processor, use 'cluster.nodes=rep("localhost",4)'. This option
%%  requires the 'rpvm' and 'snowFT' libraries. PVM must be pre-initialized
%%  and running on all nodes prior to using this option.
}
\value{
  A list of class 'chip.deconv', for which 'plot', 'print', and 'coef'
  functions exist, and containing the following elements (repeated
  'n.boot' times for each bootstrap run):
  \item{data}{The input data subset (subset selected using 'where',
    'window' and 'center'; see Description.}
  \item{fit}{Best-fit values at the locations of each probe in the input
    data.}
  \item{kernel}{Kernel used for deconvolution (provided by 'kernel'
    parameter).}
  \item{coeffs}{Non-zero coefficients (coordinate and intensity)
    for the chosen best-fit model.}
  \item{out.info}{Statistics on the best-fit solution including the LARS
    step number, BIC, RSS, etc.}
  \item{args}{All parameters input to 'chip.deconv', used for plotting
    and future reference.}
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
  deconv.entire.genome, fit.peak.profile, generate.fake.data,
  generate.binding.profile, MeDiChI-data,
  \code{lars}, \code{quadprog}, \code{Matrix}
}
\examples{
## see 'help(MeDiChI)', or...
## Run the demo yourself:

data( "halo.lowres", package="MeDiChI" )

fit <- chip.deconv( data.halo.lowres, where="Chr", fit.res=30,
               center=650000, wind=20000, max.steps=100, n.boot=10,
               kernel=kernel.halo.lowres, verbose=TRUE, boot.sample.opt="case" )

plot( fit, plot.genes=TRUE, cex=0.5, cex.lab=0.8, cex.axis=0.8 ) 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}

