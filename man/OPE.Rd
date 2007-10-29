\name{OPE}

\alias{OPE}
\alias{initOPE}
\alias{resetOPE}
\alias{adjustOPE}
\alias{predictOPE}
\alias{sampleOPE}

\title{Outer product emulator}

\description{ Initialise and fit an outer-product emlator to the
multivariate evaluations of a computer model (referred to below as a
\sQuote{simulator}).  Make predictions and sample from them.}

\usage{
initOPE(gr, kappar, Gs, Ws, NIG)
resetOPE(OPE)
adjustOPE(OPE, R, Y)
predictOPE(OPE, Rp, type = c('Student-t', 'EV'), drop = TRUE)
sampleOPE(OPE, Rp, N, drop = TRUE)
}

\arguments{
  \item{gr}{Function to compute regressors for runs}
  \item{kappar}{Function to compute residual variance for runs}
  \item{Gs}{Regression matrix for the simulator outputs}
  \item{Ws}{Variance matrix for the outputs}
  \item{NIG}{List of Normal Inverse Gamma parameters for [beta, tau]; with
 names \code{m}, \code{V}, \code{a}, and \code{d}}
  \item{OPE}{Object with class \sQuote{OPE}}
  \item{R,Rp}{Object containing the input values for the runs (see, Details,
 below)}
  \item{Y}{Matrix of outputs}
  \item{type}{Type of output returned from \code{predictOPE}}
  \item{drop}{Simplify outputs when only one run?}
  \item{N}{Number of random samples in \code{sampleOPE}}
}

\details{

For flexibility, the format of \code{R} and \code{Rp} is unrestricted.
The crucial thing is that functional calls such as \code{gr(R)},
\code{kappar(R)}, and \code{kappar(R, Rp)} all result in matrices of
the appropriate size.  So, for example, \code{R} might be a vector, a
matrix, or a dataframe, as long as \code{gr} and \code{kappar} are
written accordingly.

When predicting, the full distribution over all runs in \code{Rp} is
returned.  Often the inter-run covariances are not required, ie the
predictions are required marginally, one run at a time.  In this case
the user should call \code{predictOPE} one run at a time.  Note that
this cannot be automated within \code{predictOPE} because it would
require some knowledge of the format of \code{Rp}.

\strong{Names:} The primary source of output names are the column
names of \code{Y}, passed into \code{adjustOPE}.  Where this has not
been called (eg where \code{predictOPE} is called after
\code{initOPE}), or where \code{NULL}, output names are taken from the
row names of \code{Gs}.  Run names are taken from the row names of
\code{gr(Rp)}.

}

\value{

\item{\code{initOPE}}{Creates an object of class \sQuote{OPE}.  The other functions
take such an object and either modify it, or do something with it.}

\item{\code{resetOPE}}{Strips the ensemble from the \sQuote{OPE} object, and returns
the object.}

\item{\code{adjustOPE}}{Adds an ensemble to the \sQuote{OPE} object, and returns the
object. Only one ensemble can be added.  Multiple ensembles should be
concatenated.}

\item{\code{predictOPE}}{Predicts the simulator response at a given set of
runs.  If \code{type = "Student-t"} the prediction takes the form of a
list with three components: \code{mu}, the mean matrix; \code{Sigma}
the scale array; \code{df}, a scalar degrees of freedom.  These are
parameters of a multivariate Student-t distribution. If \code{type =
"EV"} the list summarises this information in terms of: \code{mu}, the
mean matrix; \code{Sigma}, the variance array.

Note that these outputs have their natural shape.  The mean matrix is
\eqn{n' \times q}{np by q}, where \eqn{n'}{np} is the number of rows
of \code{Rp} and \eqn{q} is the number of columns of \code{Y}.  The
variance array is four-dimensional.  The exception is when \eqn{n' =
1}{np = 1} and \code{drop = TRUE}, in which case the mean is a vector
along the outputs, and the variance is a matrix.}

\item{\code{sampleOPE}}{Samples the predictions of the simulator.  The
resulting matrix is \eqn{N \times n' \times q}{N by np by q}.  When
\eqn{n' = 1}{np = 1} and \code{drop = TRUE}, the second extent is
dropped.}

}

\references{J.C. Rougier, 2007, Efficient emulators for multivariate
deterministic functions, unpublished, available at
\url{http://www.maths.bris.ac.uk/~mazjcr/OPemulator.pdf}}

\author{Jonathan Rougier, \email{j.c.rougier@bristol.ac.uk}}

\seealso{\code{\link{rmvt}} for sampling from a multivariate Student-t.}

\examples{
#### A simple example where x in [0, 1] and theta in [0, 2pi].  This
#### example has gr() a periodic function of theta, just for fun.

## here's the true function, and we'll generate some data

g <- function(x, theta) {
  outer(x, theta, function(x, t) exp(-x) * sin(t - pi * x))
}
x <- c(0.1, 0.5, 0.7)
theta <- c(0, 1, 2, 3.5, 5) # uneven spacing more interesting
Y <- g(x, theta)

## little picture

matplot(theta, t(Y), xlim = c(0, 2 * pi), ylim = c(-1, 1),
  type = 'p', pch = 1:3, bty = 'n',
  xlab = "Theta in [0, 2pi]", ylab = "g(x, theta)",
  main = "True function and evaluations")
tfull <- seq(from = 0, to = 2*pi, len = 101)
matplot(tfull, t(g(x, tfull)), type = "l", lty = 2, add = TRUE)
legend('topright', legend = paste('x =', x),
  col = 1:3, lty = 2, pch = 1:3, bty = 'n')

## Set up the regressors and variance functions: polynomials for the
## runs regressors (should be Legendre polynomials really, shifted
## onto [0, 1]); Fourier terms for outputs regressors; power
## exponential for the runs variance function; circular correlation
## for the outputs variance matrix (note that pi cannot be too small
## or this variance is singular)

## put rownames in gr and on Gs, just for clarity

gr <- function(x) {
  robj <- cbind(1, 2*x - 1, x^2)
  rownames(robj) <- paste('x', seq(along = x), sep = '')
  robj
}

kappar <- function(x, xp = x, range = 0.5)
  exp(-abs(outer(x, xp, '-') / range)^(3/2))

Gs <- cbind(1, sin(theta), cos(theta), sin(2 * theta), cos(2 * theta))
rownames(Gs) <- paste('th', seq(along = theta), sep = '')

circular <- function(ang1, ang2 = ang1, range = pi / 1.1)
{
  smallestAngle <- function(a, b) {
    dd <- outer(a, b, '-')
    pmin(abs(dd), abs(2*pi + dd), abs(dd - 2*pi))
  }

  angles <- smallestAngle(ang1, ang2)
  ifelse(angles < range, 1 - angles / range, 0)
}

Ws <- circular(theta)

## Set up a minimal prior for the NIG (in general, thought is required
## here!)

local({
  vr <- length(gr(0))
  vs <- ncol(Gs)
  m <- rep(0, vr * vs)
  V <- diag(1^2, vr * vs)
  a <- 1
  d <- 1^2
  NIG <<- list(m = m, V = V, a = a, d = d)
})

## Now we're ready to initialise our OPE

myOPE <- initOPE(gr = gr, kappar = kappar, Gs = Gs, Ws = Ws, NIG = NIG)
 
xnew <- 0.4
pp0 <- predictOPE(myOPE, Rp = xnew) # prior prediction

## Adjust with the evaluations

myOPE <- adjustOPE(myOPE, R = x, Y = Y)

## Sanity check: predict the points we already have

pp1 <- predictOPE(myOPE, R = x)
stopifnot(
  all.equal.numeric(pp1$mu, Y, check.attributes = FALSE),
  all.equal.numeric(pp1$Sigma, array(0, dim(pp1$Sigma)))
) # phew!

## Make a prediction at some new x values, and add to the plot as
## error bars

pp2 <- predictOPE(myOPE, Rp = xnew)
pp2$mu <- c(pp2$mu) # reshape for convenience
dim(pp2$Sigma) <- rep(length(pp2$mu), 2) # 

mu <- pp2$mu
sig <- sqrt(diag(pp2$Sigma))
arrows(theta, mu + sig * qt(0.025, df = pp2$df),
       theta, mu + sig * qt(0.975, df = pp2$df),
  code = 3, angle = 90, length = 0.1, col = 'blue')
lines(tfull, g(xnew, tfull), col = 'blue')

## Add on some sampled values, interpolated using splines

rsam <- sampleOPE(myOPE, Rp = xnew, N = 10)
if (require(splines)) {
  for (i in 1:nrow(rsam)) {
    pispl <- periodicSpline(theta, rsam[i, ], period = 2*pi)
    lines(predict(pispl, tfull), col = 'darkgrey')
  }
  legend('topleft', legend = c(paste('x =', xnew, '(predicted)'), 'sampled'),
    col = c('blue', 'darkgrey'), lty = 1, pch = NA, bty = 'n')
}

## A more complicated prediction

xnew <- c(xnew, 0.8)
pp3 <- predictOPE(myOPE, Rp = xnew, type = 'EV')
}

\keyword{multivariate}
