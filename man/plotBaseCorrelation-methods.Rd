\name{plotBaseCorrelation}
\alias{plotBaseCorrelation}
\title{plotBaseCorrelation}
\usage{
plotBaseCorrelation(BaseCorrelation, ...)
}
\arguments{
  \item{BaseCorrelation}{\code{linkS4class{BaseCorrelation}}
  object}
}
\description{
Reads a BaseCorrelation object and creates a plot showing
correlation between the smooth (top right) and detail
(bottom left) wavelet coeffcients at each scale, between
each pair of variables. The relationships between
inter-pulse duration, the number of inserts, and each
individual nucleotide base are all examined.
}

