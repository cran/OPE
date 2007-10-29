\name{internals}

\alias{is.variance}
\alias{is.scalar}
\alias{is.numeric.vector}


\title{Utility functions for package OPE}

\description{Utility functions for package OPE}

\usage{
is.numeric.vector(x)
is.scalar(x)
is.variance(x, strict = TRUE, tol = 1e-06)
}

\arguments{
  \item{x}{An object of some kind}
  \item{strict}{Require \code{x} to be non-singlar?}
  \item{tol}{Numerical tolerance for spectral decomposition}
}

\value{Logical.}

\author{Jonathan Rougier \email{j.c.rougier@bristol.ac.uk}}

\examples{
is.numeric.vector(3:6) # yes
is.numeric.vector(LETTERS) # no
is.scalar(3:5) # no
is.scalar('a') # no
is.variance(matrix(c(1, 1, 1, 4), 2, 2)) # yes
is.variance(matrix(c(1, 1, 1, 1), 2, 2)) # no
is.variance(matrix(c(1, 1, 1, 1), 2, 2), strict = FALSE) # yes
}

\keyword{internal}
