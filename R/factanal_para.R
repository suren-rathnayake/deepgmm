factanal_para <- function(data, s, k, r, i, numobs) {

  lst <- list(w = list(), H = list(), mu = list(), psi = list(),
                                                   psi.inv = list())
  psi <- psi.inv <- array(0, c(k[i], r[i], r[i]))
  H <- array(0, c(k[i], r[i], r[i + 1]))
  mu <- matrix(0, r[i], k[i])

  z <- matrix(NA, nrow = numobs, ncol = r[i + 1])
  for (j in 1 : k[i]) {

    indices <- which(s == j)
    stima <- try(factanal(data[indices, ], r[i + 1], rotation = "none",
                scores = "Bartlett"), silent = TRUE)

    if (is.character(stima)) {

      psi[j,, ] <- 0.1 * diag(r[i])
      psi.inv[j,, ] <- diag(r[i])
      H[j,, ] <- matrix(runif(r[i] * r[i + 1]), r[i], r[i + 1])

      zt <- try(princomp(data[indices, ])$scores[, 1 : r[i + 1]],
                silent = TRUE)

      if (!is.character(zt)) {
        zt <- matrix(zt, ncol = r[i+1])
      } else {
        zt <- matrix(data[indices, sample(1 : r[i + 1])], ncol = r[i + 1])
      }
      z[indices, ] <- zt

    } else {

      psi[j,, ] <- diag(stima$uniq)
      H[j,, ] <- stima$load
      psi.inv[j,, ] <- diag(1/stima$uniq)
      z[indices, ] <- stima$scores
    }

    mu[, j] <- colMeans(data[indices,, drop = FALSE])
  }

  w <- matrix(table(s) / numobs)
  lst$w[i] <- list(w)
  lst$H[i] <- list(H)
  lst$mu[i] <- list(mu)
  lst$psi[i] <- list(psi)
  lst$psi.inv[i] <- list(psi.inv)

  lst
}
