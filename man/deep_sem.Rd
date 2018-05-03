\name{deep_sem}
\alias{deep_sem}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
deep_sem(y, layers, k = rep(2, layers), r = rep(1, layers), it = 50, eps = 0.001, seed = 7, init = "kmeans")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{layers}{
%%     ~~Describe \code{layers} here~~
}
  \item{k}{
%%     ~~Describe \code{k} here~~
}
  \item{r}{
%%     ~~Describe \code{r} here~~
}
  \item{it}{
%%     ~~Describe \code{it} here~~
}
  \item{eps}{
%%     ~~Describe \code{eps} here~~
}
  \item{seed}{
%%     ~~Describe \code{seed} here~~
}
  \item{init}{
%%     ~~Describe \code{init} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (y, layers, k = rep(2, layers), r = rep(1, layers), 
    it = 50, eps = 0.001, seed = 7, init = "kmeans") 
{
    ptm <- proc.time()
    a <- paste("Fit the model: seed=", seed, " k = (", paste(k, 
        collapse = " "), ") r = (", paste(r, collapse = " "), 
        ")\n", sep = "")
    cat(a)
    set.seed(seed)
    y <- scale(y)
    numobs <- nrow(y)
    p <- ncol(y)
    r <- c(p, r)
    w.list <- NULL
    H.list <- NULL
    psi.list <- NULL
    psi.list.inv <- NULL
    mu.list <- NULL
    z.list <- NULL
    for (i in 1:layers) {
        if (i == 1) {
            data <- y
        }
        else {
            data <- z[, 1:r[i], drop = FALSE]
        }
        if (init == "kmeans") {
            if (k[i] > 1) {
                s <- kmeans(data, k[i], iter.max = 100, nstart = 30, 
                  algorithm = "Hartigan-Wong")$cluster
            }
            else {
                s <- rep(1, numobs)
            }
        }
        if (init == "hclass") {
            if (k[i] > 1) {
                s <- hclass(hc(modelName = "VVV", data = y), 
                  k[i])
            }
            else {
                s <- rep(1, numobs)
            }
        }
        if (init == "random") {
            if (k[i] > 1) {
                s <- sample(1:k[i], numobs, replace = TRUE)
            }
            else {
                s <- rep(1, numobs)
            }
        }
        for (j in 1:k[i]) {
            if ((table(s)[j]) < 2) {
                s[sample(1:numobs, 2, replace = FALSE)] <- j
            }
        }
        psi <- psi.inv <- array(0, c(k[i], r[i], r[i]))
        H <- array(0, c(k[i], r[i], r[i + 1]))
        mu <- matrix(0, r[i], k[i])
        z <- NULL
        for (j in 1:k[i]) {
            stima <- try(factanal(data[s == j, ], r[i + 1], rotation = "none", 
                scores = "Bartlett"), silent = TRUE)
            if (is.character(stima)) {
                psi[j, , ] <- 0.1 * diag(r[i])
                psi.inv[j, , ] <- diag(r[i])
                H[j, , ] <- matrix(runif(r[i] * r[i + 1]), r[i], 
                  r[i + 1])
                zt <- try(princomp(data[s == j, ])$scores[, 1:r[i + 
                  1]], silent = TRUE)
                if (!is.character(zt)) {
                  zt <- matrix(zt, ncol = r[i + 1])
                }
            }
            if (is.character(zt)) {
                zt <- matrix(data[s == j, sample(1:r[i + 1])], 
                  ncol = r[i + 1])
                z <- rbind(z, zt)
            }
            if (!is.character(stima)) {
                psi[j, , ] <- diag(stima$uniq)
                H[j, , ] <- stima$load
                psi.inv[j, , ] <- diag(1/stima$uniq)
                z <- rbind(z, stima$scores)
            }
            mu[, j] <- colMeans(data[s == j, , drop = FALSE])
        }
        w <- matrix(table(s)/numobs)
        w.list[i] <- list(w)
        H.list[i] <- list(H)
        mu.list[i] <- list(mu)
        psi.list[i] <- list(psi)
        psi.list.inv[i] <- list(psi.inv)
        z.list[i] <- list(aperm(array(z[, 1:r[i + 1]], c(numobs, 
            r[i + 1], k[i])), c(3, 1, 2)))
    }
    if (layers == 1) {
        out <- deep.sem.alg.1(y, numobs, p, r[2], k, H.list, 
            psi.list, psi.list.inv, mu.list, w.list, z.list, 
            it, eps)
    }
    if (layers == 2) {
        out <- deep.sem.alg.2(y, numobs, p, r, k, H.list, psi.list, 
            psi.list.inv, mu.list, w.list, z.list, it, eps)
    }
    if (layers == 3) {
        out <- deep.sem.alg.3(y, numobs, p, r, k, H.list, psi.list, 
            psi.list.inv, mu.list, w.list, z.list, it, eps)
    }
    H <- out$H
    w <- out$w
    mu <- out$mu
    psi <- out$psi
    lik <- out$likelihood
    s <- out$s
    bic <- out$bic
    h <- out$h
    aic <- out$aic
    icl.bic <- out$icl.bic
    clc <- out$clc
    output <- list(H = H.list, w = w.list, mu = mu.list, psi = psi.list, 
        lik = lik, bic = bic, aic = aic, clc = clc, s = s, icl.bic = icl.bic, 
        h = h, k = k, r = r, numobs = numobs, elapsed.time = proc.time() - 
            ptm, seed = seed)
    message("Estimation Details:")
    cat(paste("Log-Likelihood", round(lik[length(lik)], 2), "BIC:", 
        round(bic, 2), "AIC:", round(aic, 2), "\n"))
    cat("\n")
    invisible(output)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
