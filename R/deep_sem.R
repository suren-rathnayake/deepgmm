## needed libraries

library(MASS)
library(corpcor)
library(mvtnorm)

###############################################################################
# The algorithm implements Deep Gaussian mixture models with a numer of layers
# h=1,2 and 3 layers and values for k and r.
# The case h=1 corresponds to mixtures of factor analyzers.
# See Viroli, C. & McLachlan, G.J. Stat Comput (2017).
# https://doi.org/10.1007/s11222-017-9793-z
#
# Input:
# y -  data matrix of dimension n x p (num of observations x num of variables)
# layers - the admitted values are 1, 2 or 3
#          (for the package: please insert a check that k cannot be any other
#          values and stop if for instance k>3)
# k - number of groups in the different layers.
#     a vector of integers of length layers
# r - dimensions at the different layers. By contraint it is decreasing.
#     For instance with h=3 layers and p = 10,
#     it can be r = c(9, 5, 1) but not r=c(9, 1, 5) or r=c(10, 5, 1).
#     The contraint is p > r1 > r2 ... >=1
# it - number of iterations of the EM algorithm
# eps - tolerance of the relative increment of the log-likelihod.
#       When < eps the algorithm stops.
# init - initialitation of the partition by 'kmeans' k-means (default);
#        'hclass' hierarchical clustering or 'random'
###############################################################################

deep_sem <- function(y, layers, k = rep(2, layers), r = rep(1, layers),
 it = 50, eps = 0.001, seed = 7, init = 'kmeans')
{

ptm <- proc.time()
a <- paste("Fit the model: seed=", seed, " k = (",paste(k, collapse=" "),
     ") r = (",paste(r, collapse = " "), ")\n", sep = "")
cat(a)

set.seed(seed)
y <- scale(y)
numobs <- nrow(y)
p <- ncol(y)
r <- c(p, r)

#init
w.list <- NULL
H.list <- NULL
psi.list <- NULL
psi.list.inv <- NULL
mu.list <- NULL
z.list <- NULL

for (i in 1 : layers) {

  if (i == 1) {
    data <- y
  } else {
    data <- z[, 1 : r[i], drop = FALSE]
  }

  if (init == 'kmeans') {
    if ( k[i] > 1) {
      s <- kmeans(data, k[i], iter.max = 100, nstart = 30,
              algorithm = "Hartigan-Wong")$cluster
    } else {
      s <- rep(1,numobs)
    }
  }

  if (init == 'hclass') {
    if (k[i] > 1) {
      s <- hclass(hc(modelName = "VVV", data = y),k[i])
    } else {
      s <- rep(1, numobs)
    }
  }

  if (init == 'random') {
    if (k[i] > 1) {
      s <- sample(1 : k[i], numobs, replace = TRUE)
    } else {
      s <- rep(1, numobs)
    }
  }

  for  (j in 1 : k[i]) {
    if ((table(s)[j]) < 2) {
      s[sample(1 : numobs, 2, replace = FALSE)] <- j
    }
  }

  psi <- psi.inv <- array(0, c(k[i], r[i], r[i]))
  H <- array(0, c(k[i], r[i], r[i+1]))
  mu <- matrix(0, r[i], k[i])
  z <- NULL

  for (j in 1 : k[i]) {
    stima <- try(factanal(data[s == j, ], r[i+1], rotation="none",
                 scores="Bartlett"), silent=TRUE)

    if (is.character(stima)) {
      psi[j,, ] <- 0.1 * diag(r[i])
      psi.inv[j,, ] <- diag(r[i])
      H[j,,] <- matrix(runif(r[i] * r[i+1]), r[i], r[i+1])
      zt <- try(princomp(data[s == j, ])$scores[, 1 : r[i+1]], silent = TRUE)
      if (!is.character(zt)) {
        zt <- matrix(zt, ncol = r[i+1])
      }
    }
    if (is.character(zt)) {
      zt <- matrix(data[s == j, sample(1 : r[i+1])], ncol = r[i+1])
      z <- rbind(z,zt)
    }
    if (!is.character(stima)) {
      psi[j,, ] <- diag(stima$uniq)
      H[j,, ] <- stima$load
      psi.inv[j,,] <- diag(1/stima$uniq)
      z <- rbind(z, stima$scores)
    }
    mu[, j] <- colMeans(data[s == j,, drop = FALSE])
  }

  w <- matrix(table(s)/numobs)
  w.list[i] <- list(w)
  H.list[i] <- list(H)
  mu.list[i] <- list(mu)
  psi.list[i] <- list(psi)
  psi.list.inv[i] <- list(psi.inv)
  z.list[i] <- list(aperm(array(z[, 1 : r[i+1]], c(numobs, r[i+1], k[i])),
                 c(3, 1, 2)))
 }

##############################################################################
if (layers == 1) {
  out <- deep.sem.alg.1(y, numobs, p, r[2], k, H.list, psi.list,
                        psi.list.inv, mu.list, w.list, z.list, it, eps)
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

output <- list(H = H.list, w = w.list, mu = mu.list, psi = psi.list, lik = lik,
               bic = bic, aic = aic, clc = clc, s = s, icl.bic = icl.bic,
               h = h, k = k, r = r, numobs = numobs,
               elapsed.time = proc.time() - ptm, seed = seed)

message("Estimation Details:")
cat(paste("Log-Likelihood", round(lik[length(lik)], 2), "BIC:",
      round(bic, 2), "AIC:", round(aic, 2), "\n"))
cat("\n")
invisible(output)
}
## s is the final classification
