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

  \item{shiftWindow}{}

  \item{totalTime}{}

  \item{filterNumber}{}

  \item{shrink}{}
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

