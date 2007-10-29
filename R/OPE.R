#### OPE: outer product emulator

## Code by Jonathan Rougier, j.c.rougier@bristol.ac.uk.  See
## J.C. Rougier, 2007, Efficient emulators for Multivariate
## Deterministic Funcions, under submission, currently available at
## http://www.maths.bris.ac.uk/~mazjcr/OPemulator.pdf.

is.numeric.vector <- function(x) is.vector(x) && is.numeric(x) && !any(is.na(x))

is.scalar <- function(x) is.numeric.vector(x) && length(x) == 1

is.variance <- function(x, strict = TRUE, tol = 1e-6) {
  result <- is.numeric.vector(c(x)) && is.matrix(x) && nrow(x) == ncol(x)
  result <- result &&
    all.equal.numeric(x, t(x), check.attributes = FALSE) # not bothering with tol here
  if (result) {
    ev <- eigen(x, symmetric = TRUE, only.values = TRUE)$values
    result <- all(ev / abs(ev[1]) >= (2 * strict - 1) * tol)
  }
  result
}

## Initialise with primitive objects:

## gr              function to compute regressors for the runs
## kappar          function to compute residual variance for the runs
## Gs              regression matrix for the outputs
## Ws              variance matrix for the outputs
## NIG             Normal Inverse Gamma parameters for [beta, tau], {m, V, a, d}

initOPE <- function(gr, kappar, Gs, Ws, NIG) {

  ## possible init from another OPE

  argnms <- names(formals(initOPE))
  
  if (inherits(gr, 'OPE')) {
    tmp <- gr # so we don't overwrite gr
    invisible(lapply(argnms, function(nm) assign(nm, tmp[[nm]])))
  }
  
  ## initial check of the arguments
  
  stopifnot(is.function(gr),
            is.function(kappar),
            is.matrix(Gs),
            is.variance(Ws),
            is.list(NIG),
            all(c('m', 'V', 'a', 'd') %in% names(NIG))
            )

  stopifnot(is.numeric.vector(NIG$m),
            is.variance(NIG$V),
            length(NIG$m) == nrow(NIG$V),
            is.scalar(NIG$a), # could be negative, eg OLS (yuck!)
            is.scalar(NIG$d) && NIG$d >= 0
            )

  ## package and return

  OPE <- lapply(structure(argnms, names = argnms), function(nm) get(nm))
  class(OPE) <- 'OPE'

  return(OPE)
}

resetOPE <- function(OPE) {

  if (!inherits(OPE, 'OPE'))
    stop("Expecting an object of class 'OPE'")

  gotEnsemble <- !is.null(OPE$adjust)
  if (!gotEnsemble)
    warning('Nothing to reset!')
  else
    OPE$adjust <- NULL # deleted
  return(OPE)
}

adjustOPE <- function(OPE, R, Y) {

  if (!inherits(OPE, 'OPE'))
    stop("Expecting an object of class 'OPE'")

  lapply(names(OPE), function(nm) assign(nm, OPE[[nm]])) # into the function frame
  
  ## checking the arguments

  Gr <- try(gr(R))
  if (inherits(Gr, 'try-error') || !is.matrix(Gr))
    stop("Cannot call 'gr' with argument 'R'")

  Wr <- try(kappar(R))
  if (inherits(Wr, 'try-error') || !is.variance(Wr))
    stop("Cannot call 'kappar' with argument 'R'")

  n <- nrow(Gr)
  q <- nrow(Gs)
  stopifnot(is.matrix(Y), identical(dim(Y), c(n, q)))

  vr <- ncol(Gr)
  vs <- ncol(Gs)
  v <- vr * vs
  stopifnot(length(NIG$m) == v)

  ## some obects that we will need below and later

  M <- matrix(NIG$m, vr, vs)
  C <- Y - Gr %*% tcrossprod(M, Gs)
  QWr <- chol(Wr); QWs <- chol(Ws)
  iWr <- chol2inv(QWr); iWs <- chol2inv(QWs)

  GtiWG <- kronecker(crossprod(backsolve(QWs, Gs, transpose = TRUE)),
                     crossprod(backsolve(QWr, Gr, transpose = TRUE)))
  Hr <- crossprod(Gr, iWr); Hs <- crossprod(Gs, iWs)
  GtiWc <- c(Hr %*% tcrossprod(C, Hs))

  Q <- chol(chol2inv(chol(NIG$V)) + GtiWG)
  iQall <- backsolve(Q, cbind(GtiWG, GtiWc), transpose = TRUE)
  iQGtiWG <- iQall[, 1:v, drop=FALSE]
  iQGtiWc <- iQall[, v+1]  

  mhd <- sum(C * (iWr %*% tcrossprod(C, iWs))) - drop(crossprod(iQGtiWc))
  
  ## update the NIG parameters; Eb inherits names from NIG$m if it has
  ## them

  mnew <- NIG$m + drop(NIG$V %*% (GtiWc - crossprod(iQGtiWG, iQGtiWc)))
  Vnew <- chol2inv(Q) # D^{-1}
  anew <- NIG$a + length(Y)
  dnew <- NIG$d + mhd

  ## save these, and all the stuff that will be useful in predict (in other)
  
  adjust <- list(R = R, Y = Y)
  adjust$NIG <- list(m = mnew, V = Vnew, a = anew, d = dnew)

  keep <- c('C', 'QWr', 'iWr', 'Hr', 'Q', 'iQGtiWG', 'iQGtiWc')
  adjust$other <- lapply(structure(keep, names = keep), function(nm) get(nm))

  OPE$adjust <- adjust
  return(OPE)
}

predictOPE <- function(OPE, Rp, type = c('Student-t', 'EV'), drop = TRUE) {

  if (!inherits(OPE, 'OPE'))
    stop("Expecting an object of class 'OPE'")

  lapply(names(OPE), function(nm)
         assign(nm, OPE[[nm]], pos = sys.frame(-2))) # into the function frame
  
  type <- match.arg(type)

  ## some objects that we will need below

  Grp <- try(gr(Rp))
  if (inherits(Grp, 'try-error') || !is.matrix(Grp))
    stop("Cannot call 'gr' with argument 'Rp'")

  Gp <- kronecker(Gs, Grp)

  Wrpp <- try(kappar(Rp))
  if (inherits(Wrpp, 'try-error') || !is.variance(Wrpp, strict = FALSE))
    stop("Cannot call 'kappar' with argument 'Rp'")
  Wpp <- kronecker(Ws, Wrpp)

  np <- nrow(Grp)
  q <- nrow(Gs)

  ## two possibilities: to predict from the prior or from the adjusted

  gotEnsemble <- !is.null(OPE$adjust)

  if (!gotEnsemble) {

    ## this one is easy!
    
    EfRp <- drop(Gp %*% NIG$m)
    dim(EfRp) <- c(np, q)
    
    SfRp <- crossprod(tcrossprod(chol(NIG$V), Gp)) + Wpp
    dim(SfRp) <- rep(c(np, q), 2)
    
    ap <- NIG$a
    dp <- NIG$d

  } else {
    
    ## in this case we also have to compute the mean and variance of
    ## e(Rp) and the covariance between beta and e(Rp)

    Wrp <- try(kappar(adjust$R, Rp))
    if (inherits(Wrp, 'try-error'))
      stop("Cannot call 'kappar' with arguments 'R' and 'Rp'")
  
    lapply(names(adjust$other), function(nm)
           assign(nm, adjust$other[[nm]], pos = sys.frame(-2))) # for simplicity

    iWWp <- kronecker(diag(1, q), iWr %*% Wrp)
    WpiWWp <- kronecker(Ws, crossprod(backsolve(QWr, Wrp, transpose = TRUE)))

    GtiWWp <- kronecker(t(Gs), Hr %*% Wrp)        
    iQGtiWWp <- backsolve(Q, GtiWWp, transpose = TRUE)

    Eep <- drop(crossprod(iWWp, c(C)) - crossprod(iQGtiWWp, iQGtiWc))
    Sep <- Wpp - WpiWWp + crossprod(iQGtiWWp)
    Cbep <- -NIG$V %*% (GtiWWp - crossprod(iQGtiWG, iQGtiWWp))

    ## assemble the prediction

    EfRp <- drop(Gp %*% adjust$NIG$m) + Eep
    dim(EfRp) <- c(np, q)

    tmp <- Gp %*% Cbep
    SfRp <- crossprod(tcrossprod(chol(adjust$NIG$V), Gp)) + Sep + tmp + t(tmp)
    dim(SfRp) <- rep(c(np, q), 2)
      
    ap <- adjust$NIG$a
    dp <- adjust$NIG$d

  } # end of if (!gotEnsemble)

  ## handle drop and dimnames

  nms <- rownames(Gs)
  if (gotEnsemble && !is.null(nnms <- colnames(adjust$Y)))
    nms <- nnms

  if (drop && np == 1) {

    EfRp <- EfRp[1, ]
    dim(SfRp) <- rep(length(EfRp), 2)    

    if (!is.null(nms))
      names(EfRp) <- nms
    
  } else {

    if (!is.null(nms))
      colnames(EfRp) <- nms
    if (!is.null(nms <- rownames(Grp)))
      rownames(EfRp) <- nms
  }
  
  ## repackage according to type and return; remember to check for whether moments exist

  if (type == 'Student-t')
    return(list(mu = EfRp, Sigma = (dp / ap) * SfRp, df = ap))
  else if (type == 'EV') {
    mu <- EfRp
    if (ap <= 1)
      mu[] <- NA
    Sigma <- (dp / (ap - 2)) * SfRp
    if (ap <= 2)
      Sigma[] <- NA
    return(list(mu = mu, Sigma = Sigma))
  } else stop('Never get here.')
}

sampleOPE <- function(OPE, Rp, N = 10, drop = TRUE) {

  if (!inherits(OPE, 'OPE'))
    stop("Expecting an object of class 'OPE'")
  pp <- predictOPE(OPE, Rp, type = 'Student-t', drop = FALSE)

  ## grab the dimensions before sending off to rmvt

  dd <- dim(pp$mu)
  mu <- c(pp$mu)
  Sigma <- matrix(pp$Sigma, prod(dd), prod(dd))

  rsam <- rmvt(N, mu = mu, Sigma = Sigma, df = pp$df, drop = FALSE)
  dim(rsam) <- c(N, dd)

  ## handle drop and dimnames.  Assign dimnames, and then remove if
  ## redundant

  nms <- colnames(pp$mu)
  q <- dd[2]

  if (dd[1] == 1 && drop) {
    dim(rsam) <- c(N, q)
    dimnames(rsam) <- vector("list", 2)
    if (!is.null(nms))
      colnames(rsam) <- nms
  } else {
    dimnames(rsam) <- vector("list", 3)
    if (!is.null(nms))
      dimnames(rsam)[[3]] <- nms
    if (!is.null(nms <- rownames(pp$mu)))
      dimnames(rsam)[[2]] <- nms
  }
  
  if (all(sapply(dimnames(rsam), is.null)))
    dimnames(rsam) <- NULL

  return(rsam)
}

## MV Student-t random variates

"rmvt" <- function(n = 1, mu, Sigma, df, tol = 1e-6, drop = TRUE)
{
  stopifnot(is.numeric.vector(mu),
            is.matrix(Sigma) && nrow(Sigma) == ncol(Sigma),
            length(mu) == nrow(Sigma),
            is.scalar(df))

  k <- length(mu)
  eS <- eigen(Sigma, symmetric = TRUE)
  ev <- eS$values
  if (!all(ev >= -tol * abs(ev[1])))
    stop("\'Sigma\' not non-negative definite")

  X <- matrix(rnorm(k * n), n)
  X <- tcrossprod(eS$vectors %*% diag(sqrt(pmax(ev, 0)), k), X)
  chi <- rchisq(n, df = df)
  X <- drop(mu) + sweep(X, 2, sqrt(chi/df), "/")

  nm <- names(mu)
  if (is.null(nm))
    nm <- dimnames(Sigma)[[1]]
  if (!is.null(nm))
    dimnames(X) <- list(nm, NULL)
  if (n == 1 && drop)
    drop(X)
  else t(X)
}
