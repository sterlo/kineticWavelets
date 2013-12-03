\docType{methods}
\name{wave128Window}
\alias{wave128Window}
\alias{wave128Window,KineticWavelets-method}
\title{wave128Window generic}
\usage{
wave128Window(KineticWavelets, DNAPattern, ...)

\S4method{wave128Window}{KineticWavelets}(KineticWavelets, DNAPattern,
  maxReads = 1000, shiftWindow = 64, totalTime = TRUE, filterNumber = 1,
  shrink = 1)
}
\arguments{
  \item{KineticWavelets}{\code{linkS4class{KineticWavelets}}
  object}

  \item{DNAPattern}{The DNA pattern to search for in the
  reads}

  \item{shiftWindow}{The shift window determines where in the 128 bp window the element of interest lies.  Default is 32 so that the element starts at base 33 in the window (i.e. shifted from 1 by 32).}

  \item{totalTime}{If set to TRUE, the analysis will include time spent on insert errors in the measurements.  The default totalTime==FALSE results in an analysis in which the time spent on insert errors is not included.}

  \item{filterNumber}{This allows the filter number used in the wavelet transform to be adjusted.  The default is filterNumber==1, the Haar wavelet transform.}

  \item{shrink}{The shrink option allows for only upper quantiles to be examined.  The default shrink==1 removes no data from the measurements.  If shrink is set below 1, then only the values for the upper quantile above "shrink" are used (e.g. shrink==0.05 returns the top 5 percent of the results).  Measurements are shrunk at each scale of the wavelet coefficients.}
}
\value{
A list containining the kinetics for matched sequences
}
\description{
Reads a KineticWavelets object and creates wavelets for
reads that match the DNAPattern
}
\author{
Sterling Sawaya \email{sterlingsawaya@gmail.com}
}

