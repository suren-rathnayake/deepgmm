compute.lik <- function(y, numobs, k, mu.list, H.list, psi.list, w.list) {

  value <- 1e-20
  layers <- length(k)
  p <- ncol(y)
  py <- matrix(0, numobs)
  tot.k <- prod(k)
  py.s <- matrix(0, numobs, tot.k)
  pys <- matrix(0, numobs, tot.k)

  k.comb <- apply(t(k), 2, function(x) 1 : x)
  if (is.list(k.comb)) {
    k.comb <- expand.grid(k.comb)
  }
  if (is.matrix(k.comb)) {
    k.comb <- expand.grid(split(t(k.comb), 1 : ncol(k.comb)))
  }
  if (prod(k) == 1) {
    k.comb <- matrix(k.comb, nrow =1)
  }

  for (i in 1 : tot.k)  {

    mu.tot <- mu.list[[1]][, k.comb[i, 1]]
    var.tot <- psi.list[[1]][k.comb[i, 1],, ]
    w.tot <- w.list[[1]][k.comb[i, 1]]

    if (layers > 1) {
      for (l in 2 : layers) {

        tot.H <- diag(p)
        for (m in 1 : (l - 1)) {
          tot.H <- tot.H %*% H.list[[m]][k.comb[i, m],, ]
        }

        mu.tot <- mu.tot + tot.H %*% mu.list[[l]][, k.comb[i, l]]
        var.tot <- var.tot + tot.H %*% (H.list[[l]][k.comb[i, l],, ] %*%
                   t(H.list[[l]][k.comb[i, l],, ]) +
                   psi.list[[l]][k.comb[i, l],, ]) %*% t(tot.H)
        w.tot <- w.tot * w.list[[l]][k.comb[i, l]]
      }
    }

    if (!is.positive.definite(var.tot)) {
      var.tot <- make.positive.definite(var.tot)
    }

    py.s[,i] <- dmvnorm(y, c(mu.tot), as.matrix(var.tot), log = TRUE)

    if (w.tot == 0) {
      w.tot <-  10^(-320)
    }
    pys[, i] <- log(w.tot) + py.s[, i]
  }

  cc <- 705 - max(pys)
  pys <- pys + cc
  py <- rowSums(exp(pys))

  ps.y <- exp(pys) / matrix(py, numobs, tot.k)
  ps.y <- ifelse(is.na(ps.y), 1/k, ps.y)
  py <- exp(-cc) * py


# pys_max <- apply(pys, 1, max)
# pys <- sweep(pys, 1, pys_max, '-')
# pys <- exp(pys)
# ps.y <- sweep(pys, 1, rowSums(pys), '/')
# py <- exp(pys_max) * rowSums(pys)

  s <- matrix(0, nrow = numobs, ncol = layers)
  ps.y.list <- NULL
  for (l in 1 : layers)  {

    psi.y <- matrix(0, nrow = numobs, ncol = k[l])
    for (i in 1 : k[l]) {

      index <- (k.comb[, l] == i)
      psi.y[, i] <- rowSums(ps.y[, index, drop = FALSE])
    }
    s[, l] <- apply(psi.y, 1, which.max)
    ps.y.list[[l]] <- psi.y
  }

return(list(py = py, py.s = py.s, ps.y = ps.y, k.comb = k.comb,
   s = s, ps.y.list = ps.y.list))
}

compute.est <- function(k, p, qq, ps.y, y, Ez, Ezz, mu) {
  value <- 1e-20
  H <- array(0, c(k,p,qq))
  psi <- psi.inv <- array(0, c(k,p,p))
  numobs <- dim(y)[1]
  w <- 0
  for (i in 1:k) {
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

   psi[i,,] <- (t((y - t(matrix(mu[, i], p, numobs))) *
                matrix(ps.y[, i], numobs, p)) %*%
                (matrix((y - t(matrix(mu[, i], p, numobs))) *
                matrix(ps.y[, i], numobs, p), ncol = p)) -
                t((y - t(matrix(mu[, i], p, numobs))) *
                matrix(ps.y[, i], numobs, p)) %*%
                (matrix(Ez[i,, ], nrow = numobs) *
                matrix(ps.y[, i], numobs, qq)) %*% t(H[i,, ]))

    psi[i,,] <- psi[i,,, drop = FALSE] / sum(ps.y[, i])
    psi[i,,] <- ifelse((psi[i,,] == 0), value,psi[i,, ])
    if (p > 1) {
      psi[i,, ] <- diag(diag(psi[i,, ]))
      psi.inv[i,, ] <- diag(1/diag(psi[i,, ]))
    }
    psi[i,, ]  <- ifelse(is.na(psi[i,, ]), 1, psi[i,, ])
    psi.inv[i,,] <- ifelse(is.na(psi.inv[i,, ]), 1, psi.inv[i,, ])
    mu[, i] <- colSums(matrix(ps.y[, i], numobs, p) *
               (y - t(matrix(H[i,, ], ncol = qq) %*%
               t(Ez[i,, ]))))/sum(ps.y[, i])
    w[i] <- mean(ps.y[, i])
  }
  return(list(H = H,mu = mu, psi = psi, psi.inv = psi.inv, w = w))
}

makeSymm <- function(m) {
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  return(m)
}

ma <- function(x, n = 5) {
  filter(x, rep(1/n, n), sides = 1)
}

entr <- function(z) {
  numobs <- nrow(z)
  numg <- ncol(z)
  temp <- 0
  z <- ifelse(z==0, z + 0.000000000000000000000001, z)
  for (i in 1 : numg)  {
    for (j in 1:numobs) {
      temp <- temp + (z[j, i] *log(z[j, i]))
    }
  }
  return(-temp)
}
