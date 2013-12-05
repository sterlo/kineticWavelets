\name{grepSeq}
\alias{grepSeq}
\title{Searches sequences for a given pattern.}
\usage{
grepSeq(seqz, patternz, returnCounts = T)
}
\arguments{
  \item{seqz}{Character vector representing the to search
  for a pattern.}

  \item{patternz}{Regular expression to search for}

  \item{returnCounts}{if \code{TRUE} print the number of
  matches.}
}
\value{
A vector containing the matches.
}
\description{
The function returns a vector the length of the original
sequence (seqz) that contains a 1 if the position matched
the pattern otherwise 0.
}
\examples{
grepSeq('ATCTCTCTTCTCTTC','TC')
# 6 regions found
# [1] 0 1 1 1 1 1 1 0 1 1 1 1 0 1 1
}
\author{
Sterling Sawaya \email{sterlingsawaya@gmail.com}
}

