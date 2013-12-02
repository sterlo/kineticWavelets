\docType{methods}
\name{baseCorrelation}
\alias{baseCorrelation}
\alias{baseCorrelation,KineticWavelets-method}
\title{baseCorrelation generic}
\usage{
baseCorrelation(KineticWavelets, ...)

\S4method{baseCorrelation}{KineticWavelets}(KineticWavelets,
  DNAPattern = NULL, minReadLength = 200, maxReads = 10000)
}
\arguments{
  \item{KineticWavelets}{\code{linkS4class{KineticWavelets}}
  object}

  \item{DNAPattern}{The DNA pattern to search for in the
  reads.}

  \item{minReadLength}{minimum read length}

  \item{maxReads}{maximum number of reads to correlate.}
}
\value{
A list containining the base correlations

A list containing all the base correlations.
}
\description{
Reads a KineticWavelets object and perform base correlation
}
\author{
Sterling Sawaya \email{sterlingsawaya@gmail.com}
}

