ppca_para <- function(data, s, k, r, i, numobs)  {

  lst <- list(w = list(), H = list(), mu = list(), psi = list(),
                                                   psi.inv = list())
  psi <- psi.inv <- array(0, c(k[i], r[i], r[i]))
  H <- array(0, c(k[i], r[i], r[i + 1]))
  mu <- matrix(0, r[i], k[i])

  z <- matrix(NA, nrow = numobs, ncol = r[i + 1])
  for (j in 1 : k[i]) {

    q <- r[i + 1]
    indices <- which(s == j)
    mu[, j] <- colMeans(data[indices,, drop = FALSE])
    Si <- cov(data[indices, ])
    psi[j,, ] <-  diag(diag(Si))
    Di.sqrt <- diag(sqrt(diag(diag(diag(Si)))))
    inv.Di.sqrt <- diag(1 / diag(Di.sqrt))
    eig.list <- eigen(inv.Di.sqrt %*% Si %*% inv.Di.sqrt)
    eigH <- eig.list$vectors
    sort.lambda <- sort(eig.list$values, decreasing = TRUE,
                                         index.return = TRUE)
    lambda <- sort.lambda$x
    ix.lambda   <- sort.lambda$ix
    sigma2 <- mean(lambda[(q + 1) : ncol(data)])
    if (q == 1) {
      H[j,, ] <- Di.sqrt %*% eigH[, ix.lambda[1 : q]] %*%
                      diag((lambda[1 : q] - sigma2), q)
      z[indices, ] <-  sweep(data[indices,, drop = FALSE], 2,
                         mu[, j, drop = FALSE], '-') %*%
                          t(1 / (t(H[j,, ]) %*% H[j,, ] +
                            diag(sigma2, q)) %*% t(H[j,, ]))
    } else {
      H[j,, ] <- Di.sqrt %*% eigH[, ix.lambda[1 : q]] %*%
                      diag((lambda[1 : q] - sigma2))

      z[indices, ] <-  sweep(data[indices,, drop = FALSE], 2,
                           mu[, j, drop = FALSE], '-') %*%
                           t(chol.inv(t(H[j,, ]) %*% H[j,, ] +
                            diag(sigma2, q)) %*% t(H[j,, ]))
    }
  }

  w <- matrix(table(s) / numobs)
  
  return(list(w = w, H = H, mu = mu, psi = psi, psi.inv = psi.inv, z = z))
}
