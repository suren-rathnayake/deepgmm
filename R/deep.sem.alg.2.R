# internal function - it performs the Stochastic EM algorithm for fitting
# the model with h=2 layers

deep.sem.alg.2 <- function(y, numobs, p, r, k, H.list, psi.list, psi.list.inv,
                           mu.list, w.list, z.list, it, eps) {
likelihood <- NULL
hh <- 0
ratio <- 1000
layers <- length(k)
#################################
#### compute the likelihood #####
#################################
out <- compute.lik(y,numobs,k,mu.list,H.list,psi.list,w.list)
py <- out$py
ps.y <- out$ps.y
ps.y.list <- out$ps.y.list
k.comb <- out$k.comb
s <- out$s
tot.k <- prod(k)
temp <- sum(log(py))
likelihood <- c(likelihood, temp)
#####################################################
while ((hh < it) & (ratio > eps )) {
  hh <- hh + 1
###############################################################
###################  first layer ##############################
###############################################################
l <- 1
yy <- y
z2 <- z2.one <- array(0, c(numobs, r[l+1], k[l], k[l+1]))
zz2 <- array(0, c(numobs, r[l + 1], r[l + 1], k[l], k[l + 1]))
z <- z.one <- array(0, c(numobs, r[l + 1], k[l]))
zz <- array(0, c(numobs, r[l + 1], r[l + 1], k[l]))

for (p1 in 1 : k[l]) {
  for (p2 in 1 : k[l + 1]) {
    A <- ginv(H.list[[l + 1]][p2,, ] %*% t(H.list[[l + 1]][p2,, ]) +
          psi.list[[l + 1]][p2,, ]) + t(H.list[[l]][p1,, ]) %*%
          (psi.list.inv[[l]][p1,, ]) %*% (H.list[[l]][p1,, ])

    b <- ginv(H.list[[l + 1]][p2,, ] %*% t(H.list[[l + 1]][p2,, ]) +
          psi.list[[l + 1]][p2,, ]) %*%
          matrix(mu.list[[l + 1]][, p2], r[l + 1], numobs) +
          t(H.list[[l]][p1,, ]) %*% (psi.list.inv[[l]][p1,, ]) %*%
          (t(yy) - matrix(mu.list[[l]][, p1], r[l], numobs))

    chsi <- ginv(A)
    if (!isSymmetric(chsi)) {
      chsi <- makeSymm(chsi)
    }
    roy <- chsi %*% b
    roy.quadro <- array(apply(roy, 2, function(x) x %*% t(x)),
                       c(r[l + 1], r[l + 1], numobs))
    zz2[,,, p1, p2] <- aperm(array(chsi, c(r[l + 1], r[l + 1], numobs)) +
                       roy.quadro, c(3, 1, 2))
    z2.one[,, p1, p2] <- rmvnorm(numobs, rep(0, r[l + 1]), chsi) + t(roy)
    z2[,,p1,p2] <- t(roy)
  }
}

for (i in 1 : k[l + 1]) {
  prob <- ps.y.list[[l + 1]][, i, drop = FALSE]
  z <- z + array(z2[,,, i, drop = FALSE] * array(rowSums(prob),
         c(numobs, r[l + 1], k[l], 1)), c(numobs, r[l + 1], k[l]))
  z.one <- z.one + array(z2.one[,,, i, drop = FALSE] *
                   array(rowSums(prob), c(numobs, r[l + 1], k[l], 1)),
                   c(numobs, r[l+1], k[l]))
  zz <- zz + array(zz2[,,,, i, drop = FALSE] *
             array(rowSums(prob), c(numobs, r[l + 1], r[l + 1], k[l], 1)),
            c(numobs, r[l + 1], r[l + 1], k[l]))
}

z.list[[l]] <- aperm(z.one, c(3, 1, 2))
out <- compute.est(k[l], r[l], r[l+1], ps.y.list[[l]], yy,
        aperm(z, c(3, 1, 2)), aperm(zz, c(4, 2, 3, 1)), mu.list[[l]])

H.list[[l]] <- out$H
psi.list[[l]] <- out$psi
psi.list.inv[[l]] <- out$psi.inv
mu.list[[l]] <- out$mu
w.list[[l]] <- out$w

###############################################################
################### second layer ##############################
###############################################################
l <- 2
yy <- matrix(0, numobs, r[l])
zz <- z.list[[l - 1]]
for (i in 1 : k[l - 1]) {
  yy <- yy + matrix(zz[i, ,, drop = FALSE] *
           array(rowSums(ps.y.list[[l - 1]][, i, drop = FALSE]),
           c(1, numobs, r[l])), numobs, r[l])
}

z <- z.one <- array(0, c(numobs, r[l + 1], k[l]))
zz <- array(0, c(numobs, r[l + 1], r[l + 1], k[l]))

for (p1 in 1 : k[l]) {
  A <- diag(r[l + 1]) + t(H.list[[l]][p1,, ]) %*%
       (psi.list.inv[[l]][p1,, ]) %*% (H.list[[l]][p1,, ])
  b <- t(H.list[[l]][p1,, ]) %*% (psi.list.inv[[l]][p1,, ]) %*%
       (t(yy) - matrix(mu.list[[l]][, p1], r[l], numobs))

  chsi <- ginv(A)

  if (!isSymmetric(chsi)) {
    chsi <- makeSymm(chsi)
  }

  roy <- chsi%*%b
  roy.quadro <- array(apply(roy, 2, function(x) x %*% t(x)),
                  c(r[l + 1], r[l + 1], numobs))
  zz[,,,p1] <- aperm(array(chsi, c(r[l + 1], r[l + 1], numobs)) +
                 roy.quadro,c(3, 1, 2))
  z.one[,, p1] <- t(roy) + rmvnorm(numobs, rep(0, r[l + 1]), chsi)
  z[,, p1] <- t(roy)
}

z.list[[l]] <- aperm(z.one, c(3, 1, 2))
out <- compute.est(k[l], r[l], r[l + 1], ps.y.list[[l]], yy,
        aperm(z, c(3, 1, 2)), aperm(zz, c(4, 2, 3, 1)), mu.list[[l]])

H.list[[l]] <- out$H
psi.list[[l]] <- out$psi
psi.list.inv[[l]] <- out$psi.inv
mu.list[[l]] <- out$mu
w.list[[l]] <- out$w

out <- compute.lik(y, numobs, k, mu.list, H.list, psi.list, w.list)
py <- out$py
ps.y <- out$ps.y
ps.y.list <- out$ps.y.list
k.comb <- out$k.comb
s <- out$s

lik <- sum(log(py))
likelihood <- c(likelihood, lik)

if (hh < 5) {
  ratio <- 2 * eps
}

if (hh > 5) {
  ratio <- (ma(likelihood, 5) [hh + 1] - ma(likelihood, 5) [hh]) /
            abs(ma(likelihood,5)[hh])
}

}

h <- 0
for (j in 1 : layers) {
  h <- h + (k[j] - 1) + (r[j] * r[j + 1]) * k[j] + r[j] * k[j] + k[j] * r[j]  #- (r[j+1]*(r[j+1]-1)/2)*k[j]
}
if (layers > 1) {
  for (j in 2:layers) {
    h <- h - (r[j] * k[j] * (r[j] - 1) / 2)
  }
}

bic <- -2*lik+h*log(numobs)
aic <- -2*lik+2*h
EN <- entr(ps.y.list[[1]])
clc <- -2*lik+2*EN
icl.bic <- -2 * lik + 2 * EN + h * log(numobs)

out <- list(H = H.list, w = w.list, mu = mu.list, psi = psi.list,
            likelihood = likelihood, bic = bic, aic = aic, clc = clc,
            s = s, icl.bic = icl.bic, h = h, ps.y = ps.y)
return(out)
}

compute.lik <- function(y, numobs, k, mu.list, H.list, psi.list, w.list) {
  value <- 1e-20
  layers <- length(k)
  p <- ncol(y)
  py <- matrix(0,numobs)
  tot.k <- prod(k)
  py.s <- matrix(0,numobs,tot.k)
  pys <- matrix(0,numobs,tot.k)
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
    mu.tot <- mu.list[[1]][,k.comb[i, 1]]
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
        w.tot <- w.tot*w.list[[l]][k.comb[i,l]]
      }
    }

    if (!is.positive.definite(var.tot)) {
      var.tot <- make.positive.definite(var.tot)
    }

    py.s[,i] <- dmvnorm(y, c(mu.tot), as.matrix(var.tot), log = TRUE)

    if (w.tot == 0) {
      w.tot <- 10^(-320)
    }
    pys[, i] <- log(w.tot) + py.s[, i]
  }

  cc <- 705 - max(pys)
  pys <- pys + cc
  py <- rowSums(exp(pys))

  ps.y <- exp(pys) / matrix(py, numobs, tot.k)
  ps.y <- ifelse(is.na(ps.y), 1/k, ps.y)
  py <- exp(-cc)*py

  s <- matrix(0, numobs, layers)
  ps.y.list <- NULL
  for (l in 1 : layers)  {
    psi.y <- matrix(0, numobs, k[l])
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
  H <- array(0,c(k,p,qq))
  psi <- psi.inv <- array(0,c(k,p,p))
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
