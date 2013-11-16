\name{MeDiChI-data}
\alias{MeDiChI-data}
\alias{data.halo.hires}
\alias{data.halo.lowres}
\alias{data.yeast.gcn4}
\alias{kernel.halo.hires}
\alias{kernel.halo.lowres}
\alias{kernel.yeast.gcn4}
\alias{kernel.yeast.gcn4.jbd}
\alias{gene.coords}
\docType{data}
\title{ChIP-chip data supporting the MeDiChI package}
\description{
  Three different ChIP-chip datasets of varying tiling resolution,
  referenced in Reference (1) below.
}
\usage{
data(halo.hires)
data(halo.lowres)
data(yeast.gcn4)
}
\format{
    The format of "halo.lowres" is:
      chr [1:3] "data.halo.lowres" "kernel.halo.lowres" "gene.coords"
      The format of "halo.hires" is:
      chr [1:3] "data.halo.hires" "kernel.halo.hires" "gene.coords"
      The format of "yeast.gcn4" is:
      chr [1:4] "data.yeast.gcn4" "kernel.yeast.gcn4"
     "kernel.yeast.gcn4.jbd" "gene.coords"
}  
\details{
 "halo.lowres" loads:\cr
    "data.halo.lowres" -- low-resolution HaloSpan Halobacterium sp. TfbD
      ChIP-chip data, formatted for input to chip.deconv() and
      deconv.entire.genome(). Data is described in ref. (1) below.\cr
    "kernel.halo.lowres" -- deconvolution kernel learned for this data
    via fit.peak.profile().\cr
    "gene.coords" -- map of gene coordinates for this organism for use
     in plotting functions.\cr 
  "halo.hires" loads:\cr
    "data.halo.hires" -- high-resolution Nimblegen format Halobacterium sp. TfbD
      ChIP-chip data, formatted for input to chip.deconv() and
      deconv.entire.genome(). Data is described in ref. (1) below.\cr
    "kernel.halo.hires" -- deconvolution kernel for this data via fit.peak.profile().\cr
    "gene.coords" -- map of gene coordinates for this organism for use
      in plotting functions.\cr
  "yeast.gcn4" loads:\cr
    "data.yeast.gcn4" -- 266bp-resolution S. cerevisiae Gcn4
      ChIP-chip data, formatted for input to chip.deconv() and
      deconv.entire.genome(). Data is described in ref. (1) below.\cr
    "kernel.yeast.gcn4" -- deconvolution kernel for this data via fit.peak.profile().\cr
    "kernel.yeast.gcn4.jbd" -- deconvolution kernel used by ref. (2) for
      the JBD algorithm.\cr
    "gene.coords" -- map of gene coordinates for this organism for use
      in plotting functions.\cr
}
\source{
  See references below.
}
\references{
  (1). Reiss, DJ and Facciotti, MT and Baliga, NS. (2007). "Model-based
  deconvolution of genome-wide DNA binding",
  Bioinformatics; doi: 10.1093/bioinformatics/btm592. \url{http://baliga.systemsbiology.net/medichi}\cr
  (2) Qi, Y and et al. (2006). "High-resolution computational models of genome
  binding events", Nature Biotechnol, 24(8), 963-970. \url{http://cgs.csail.mit.edu/jbd}.
}
\examples{
data(halo.hires)

## Then try print(data.halo.lowres) to see the format.

data(halo.lowres)
data(yeast.gcn4)
}
\keyword{datasets}
\seealso{
  chip.deconv, deconv.entire.genome, fit.peak.profile, generate.fake.data,
  generate.binding.profile
}
