\name{rmvt}

\alias{rmvt}

\title{Multivariate Student-t quantities}

\description{Generate random samples from a multivariate Student-t
with given mean vector, scale matrix, and degrees of freedom.}

\usage{
rmvt(n = 1, mu, Sigma, df, tol = 1e-06, drop = TRUE)
}

\arguments{
  \item{n}{Number of samples}
  \item{mu}{Mean vector}
  \item{Sigma}{Scale matrix}
  \item{df}{Degrees of freedom}
  \item{tol}{Numerical tolerance for spectral decomposition}
  \item{drop}{Simplify when \code{n = 1}?}
}

\value{A matrix with \code{n} rows and \code{length(mu)} columns.  A
vector if \code{n == 1 && drop}.}

\author{Jonathan Rougier \email{j.c.rougier@bristol.ac.uk}}

\examples{
mu <- c(1, 2)
Sigma <- matrix(c(1, 1, 1, 4), 2, 2)
X <- rmvt(n = 10, mu = mu, Sigma = Sigma, df = 2)
}

\keyword{distribution}
\keyword{multivariate}
