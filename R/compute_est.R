compute.est <- function(k, p, qq, ps.y, y, Ez, Ezz, mu) {

value <- .Machine$double.eps #1e-20
H <- array(0, c(k, p, qq))
psi <- psi.inv <- array(0, c(k, p, p))
numobs <- dim(y)[1]
w <- 0
for (i in 1 : k) {

  if (qq > 1) {
    EEzz <- apply((Ezz[i,,, ] * (aperm(array(t(t(ps.y[, i])),
               c(numobs, qq, qq)), c(2, 3, 1)))), 1, rowSums) / sum(ps.y[, i])
  }
  if (qq == 1) {
    EEzz <- Ezz[i,,, ] %*% (ps.y[, i]) / sum(ps.y[, i])
  }

  H[i,, ] <- (t((y - t(matrix(mu[, i], p, numobs))) *
               matrix(ps.y[,i], numobs, p)) %*%
               (matrix(Ez[i,, ], nrow = numobs) *
               matrix(ps.y[, i], numobs, qq))) %*% ginv(EEzz) / sum(ps.y[, i])

  H[i,, ] <- ifelse(is.na(H[i,, ]), mean(H[i,, ], na.rm = TRUE), H[i,, ])

  psi[i,, ] <- (t((y - t(matrix(mu[, i], p, numobs))) *
                matrix(ps.y[, i], numobs, p)) %*%
                (matrix((y - t(matrix(mu[, i], p, numobs))) *
                matrix(ps.y[, i], numobs, p), ncol = p)) -
                t((y - t(matrix(mu[, i], p, numobs))) *
                matrix(ps.y[, i], numobs, p)) %*%
                (matrix(Ez[i,, ], nrow = numobs) *
                matrix(ps.y[, i], numobs, qq)) %*% t(H[i,, ]))

  psi[i,, ] <- psi[i,,, drop = FALSE] / sum(ps.y[, i])
  psi[i,, ] <- ifelse((psi[i,, ] == 0), value, psi[i,, ])
  if (p > 1) {
    psi[i,, ] <- diag(diag(psi[i,, ]))
    psi.inv[i,, ] <- diag(1/diag(psi[i,, ]))
  }
  psi[i,, ]  <- ifelse(is.na(psi[i,, ]), 1, psi[i,, ])
  psi.inv[i,, ] <- ifelse(is.na(psi.inv[i,, ]), 1, psi.inv[i,, ])
  mu[, i] <- colSums(matrix(ps.y[, i], numobs, p) *
               (y - t(matrix(H[i,, ], ncol = qq) %*%
               t(Ez[i,, ]))))/sum(ps.y[, i])
  w[i] <- mean(ps.y[, i])
}

return(list(H = H, mu = mu, psi = psi, psi.inv = psi.inv, w = w))
}
