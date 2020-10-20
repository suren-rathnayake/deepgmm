factanal_para <- function(data, s, k, r, i, numobs) {

  lst <- list(w = list(), H = list(), mu = list(), psi = list(),
                                                   psi.inv = list())
  psi <- psi.inv <- array(0, c(k[i], r[i], r[i]))
  H <- array(0, c(k[i], r[i], r[i + 1]))
  mu <- matrix(0, r[i], k[i])

  z <- matrix(NA, nrow = numobs, ncol = r[i + 1])
  for (j in 1 : k[i]) {

    indices <- which(s == j)
    stima <- try(factanal(data[indices, ], r[i + 1], rotation = "varimax",
                scores = "regression"), silent = TRUE)

    if (is.character(stima)) {

      zt <- try(princomp(data[indices, ])$scores[, 1 : r[i + 1]],
                silent = TRUE)

      if (!is.character(zt)) {
        zt <- matrix(zt, ncol = r[i+1])
      } else {
        zt <- matrix(data[indices, sample(1 : r[i + 1])], ncol = r[i + 1])
      }
      
      H[j,, ] <- t(data[s == j, ]) %*% zt %*% ginv(t(zt) %*% zt)
      uu <- data[indices, ] - t(H[j,, ] %*% t(zt))
      psi[j,, ] <- diag(diag(var(uu)))
      psi.inv[j,, ] <-  ginv(psi[j,, ])

      z[indices, ] <- zt

    } else {

      H[j,, ] <- t(data[indices, ]) %*% stima$scores %*% 
                           ginv(t(stima$scores) %*% stima$scores)
      uu <- data[indices, ] - t(H[j,, ] %*% t(stima$scores))
      psi[j,, ] <- diag(diag(var(uu)))
      psi.inv[j,, ] <- ginv(psi[j,, ])
      z[indices, ] <- stima$scores

    }

    mu[, j] <- colMeans(data[indices,, drop = FALSE])
  }

  w <- matrix(table(s) / numobs)
  return(list(w = w, H = H, mu = mu, psi = psi, psi.inv = psi.inv, z = z))
}
