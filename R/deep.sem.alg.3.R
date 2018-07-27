deep.sem.alg.3 <- function(y, numobs, p, r, k, H.list, psi.list, psi.list.inv,
                            mu.list, w.list, it, eps) {

likelihood <- NULL
hh <- 0
ratio <- 1000
layers <- length(k)
#################################
#### compute the likelihood #####
#################################
out <- compute.lik(y, numobs, k, mu.list, H.list, psi.list, w.list)
py <- out$py
ps.y <- out$ps.y
ps.y.list <- out$ps.y.list
k.comb <- out$k.comb
s <- out$s
tot.k <- prod(k)
temp <- sum(log(py))
likelihood <- c(likelihood, temp)
#####################################################
z.list <- NULL
while ((hh < it) & (ratio > eps )) {
  hh <- hh+1
  ###############################################################
  ###################  first layer ##############################
  ###############################################################
  l <- 1
  yy <- y
  z2 <- z2.one <- array(0, c(numobs, r[l + 1], k[l], k[l + 1], k[l + 2]))
  zz2 <- array(0, c(numobs, r[l +1], r[l+1], k[l], k[l+1], k[l + 2]))
  z <- z.one <- array(0, c(numobs, r[l+1], k[l]))
  zz <- array(0, c(numobs, r[l+1], r[l+1], k[l]))

  for (p1 in 1 : k[l]) {
    for (p2 in 1 : k[l+1]) {
      for (p3 in 1 : k[l+2]) {
        sigma.tilde.inv <- ginv(H.list[[l+1]][p2,, ] %*%
                            (H.list[[l+2]][p3,,] %*% t(H.list[[l+2]][p3,,]) +
                            psi.list[[l+2]][p3,,]) %*% t(H.list[[l+1]][p2,,]) +
                            psi.list[[l+1]][p2,,])

        A <- sigma.tilde.inv + t(H.list[[l]][p1,, ]) %*%
              (psi.list.inv[[l]][p1,, ]) %*% (H.list[[l]][p1,, ])

        mu.tilde  <- matrix(mu.list[[l + 1]][, p2] + H.list[[l + 1]][p2,, ] %*%
                     (mu.list[[l + 2]][, p3]), r[l + 1], numobs)

        b <- sigma.tilde.inv %*% mu.tilde +
             t(H.list[[l]][p1,, ]) %*% (psi.list.inv[[l]][p1,, ]) %*%
             (t(yy) - matrix(mu.list[[l]][, p1], r[l], numobs))

        chsi  <- ginv(A)
        if (!isSymmetric(chsi)) {
          chsi <- makeSymm(chsi)
        }

        roy <- chsi %*% b
        roy.quadro <- array(apply(roy, 2, function(x) x %*% t(x)),
                            c(r[l + 1], r[l + 1], numobs))
        zz2[,,, p1, p2, p3]  <- aperm(array(chsi,
                                     c(r[l + 1], r[l + 1], numobs)) +
                                     roy.quadro, c(3, 1, 2))
        z2.one[,, p1, p2, p3]  <- rmvnorm(numobs, rep(0, r[l + 1]), chsi) +
                                      t(roy)
        z2[,, p1, p2, p3] <- t(roy)
      }
    }
  }

  for (i1 in 1 : k[l + 1]) {
    for (i2 in 1:k[l+2]) {
      prob <- ps.y.list[[l + 1]][, i1, drop = FALSE] *
                ps.y.list[[l + 2]][, i2, drop = FALSE]
      z <- z + array(z2[,,, i1, i2, drop = FALSE] *
                array(rowSums(prob), c(numobs, r[l + 1], k[l], 1, 1)),
                         c(numobs, r[l + 1], k[l]))
      z.one  <- z.one + array(z2.one[,,, i1, i2, drop = FALSE] *
                  array(rowSums(prob), c(numobs, r[l + 1], k[l], 1, 1)),
                  c(numobs, r[l + 1], k[l]))
      zz <- zz + array(zz2[,,,, i1, i2, drop = FALSE] *
              array(rowSums(prob), c(numobs, r[l + 1],
                r[l + 1], k[l], 1, 1)), c(numobs, r[l + 1], r[l + 1], k[l]))
    }
  }

  z.list[[l]] <- aperm(z.one, c(3, 1, 2))
  out <- compute.est(k[l], r[l], r[l + 1], ps.y.list[[l]], yy,
            aperm(z, c(3, 1, 2)), aperm(zz, c(4, 2, 3, 1)), mu.list[[l]])

  H.list[[l]] <- out$H
  psi.list[[l]] <- out$psi
  psi.list.inv[[l]] <- out$psi.inv
  mu.list[[l]] <- out$mu
  w.list[[l]] <- out$w
  ###############################################################
  ###################  second layer #############################
  ###############################################################

  l <- 2
  yy <- matrix(0, numobs, r[l])
  zz <- z.list[[l - 1]]
  for (i in 1 : k[l - 1]) {
    yy <- yy + matrix(zz[i,,, drop = FALSE] *
            array(rowSums(ps.y.list[[l - 1]][, i, drop = FALSE]),
            c(1, numobs, r[l])), numobs, r[l])
  }
  z2 <- z2.one <- array(0, c(numobs, r[l + 1], k[l], k[l + 1]))
  zz2 <- array(0, c(numobs, r[l + 1], r[l + 1], k[l], k[l + 1]))
  z <- z.one <- array(0, c(numobs, r[l + 1], k[l]))
  zz <- array(0, c(numobs, r[l + 1], r[l + 1], k[l]))

  for (p1 in 1 : k[l]) {
    for (p2 in 1:k[l+1])  {
      A <- ginv(H.list[[l + 1]][p2,, ] %*% t(H.list[[l + 1]][p2,, ]) +
           psi.list[[l + 1]][p2,, ]) + t(H.list[[l]][p1,, ]) %*%
           (psi.list.inv[[l]][p1,, ]) %*% (H.list[[l]][p1,, ])

      b <- ginv(H.list[[l + 1]][p2,, ] %*% t(H.list[[l + 1]][p2,, ]) +
           psi.list[[l + 1]][p2,, ]) %*%
           matrix(mu.list[[l + 1]][, p2], r[l + 1], numobs) +
           t(H.list[[l]][p1,, ]) %*% (psi.list.inv[[l]][p1,, ]) %*%
           (t(yy) - matrix(mu.list[[l]][, p1], r[l], numobs))

      chsi  <- ginv(A)

      if (!isSymmetric(chsi)) {
        chsi  <- makeSymm(chsi)
      }
      roy  <- chsi %*% b
      roy.quadro  <- array(apply(roy, 2, function(x) x %*%
                          t(x)), c(r[l + 1], r[l + 1], numobs))
      zz2[,,, p1, p2]  <- aperm(array(chsi, c(r[l + 1], r[l + 1], numobs)) +
                            roy.quadro, c(3, 1, 2))
      z2.one[,, p1, p2]  <- rmvnorm(numobs, rep(0, r[l + 1]), chsi) + t(roy)
      z2[,, p1, p2]  <- t(roy)
    }
  }

  for (i in 1 : k[l + 1]) {
    prob  <- ps.y.list[[l + 1]][, i, drop = FALSE]
    z  <- z + array(z2[,,, i, drop = FALSE] *
            array(rowSums(prob), c(numobs, r[l + 1], k[l], 1)),
            c(numobs, r[l + 1], k[l]))
    z.one  <- z.one + array(z2.one[,,, i, drop = FALSE] *
              array(rowSums(prob), c(numobs, r[l + 1], k[l], 1)),
              c(numobs, r[l + 1], k[l]))
    zz <- zz + array(zz2[,,,, i, drop = FALSE] * array(rowSums(prob),
            c(numobs, r[l + 1], r[l + 1], k[l], 1)),
            c(numobs, r[l + 1], r[l + 1], k[l]))
  }

  z.list[[l]] <- aperm(z.one, c(3, 1, 2))
  out  <- compute.est(k[l], r[l], r[l +1 ], ps.y.list[[l]], yy,
            aperm(z, c(3, 1, 2)), aperm(zz, c(4, 2, 3, 1)), mu.list[[l]])

  H.list[[l]] <- out$H
  psi.list[[l]] <- out$psi
  psi.list.inv[[l]] <- out$psi.inv
  mu.list[[l]] <- out$mu
  w.list[[l]] <- out$w

  ###############################################################
  ###################  third layer ##############################
  ###############################################################
  l <- 3
  yy <- matrix(0, numobs, r[l])
  zz <- z.list[[l-1]]
  for (i in 1 : k[l - 1]) {
    yy <- yy + matrix(zz[i,,, drop = FALSE] *
          array(rowSums(ps.y.list[[l-1]][, i, drop = FALSE]),
           c(1, numobs, r[l])), numobs, r[l])
  }

  z <- z.one<-array(0, c(numobs, r[l + 1], k[l]))
  zz <- array(0,c(numobs, r[l + 1], r[l + 1], k[l]))

  for (p1 in 1:k[l]) {

   A <- diag(r[l + 1]) + t(H.list[[l]][p1,, ]) %*%
          (psi.list.inv[[l]][p1,, ]) %*% (H.list[[l]][p1,, ])

   b <- t(H.list[[l]][p1,, ]) %*% (psi.list.inv[[l]][p1,, ]) %*%
        (t(yy) - matrix(mu.list[[l]][, p1], r[l], numobs))

   chsi <- ginv(A)
   if (!isSymmetric(chsi)) {
     chsi <- makeSymm(chsi)
   }

   roy <- chsi %*% b
   roy.quadro <- array(apply(roy, 2, function(x) x %*%
                  t(x)), c(r[l + 1], r[l + 1], numobs))
   zz[,,, p1]  <- aperm(array(chsi, c(r[l + 1], r[l + 1], numobs)) +
                    roy.quadro, c(3, 1, 2))
   z.one[,, p1] <- t(roy) + rmvnorm(numobs, rep(0, r[l + 1]), chsi)
   z[,, p1]  <- t(roy)
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

print(lik)

  if (hh < 5) {
    ratio <- 2 * eps
  }
  if (hh > 5) {
    ratio <- (ma(likelihood, 5)[hh + 1] - ma(likelihood, 5)[hh]) /
              abs(ma(likelihood, 5)[hh])
  }
}

h <- 0
for (j in 1 : layers) {
  h <- h + (k[j] - 1) + (r[j] * r[j + 1]) * k[j] + r[j] * k[j] + k[j] * r[j]
}

if (layers > 1) {
  for (j in 2 : layers) {
    h <- h - (r[j] * k[j] * (r[j] - 1)/2)
  }
}

bic <- -2 * lik + h * log(numobs)
aic <- -2 * lik + 2 * h
EN <- entr(ps.y.list[[1]])
clc <- -2 * lik + 2 * EN
icl.bic <- -2 * lik + 2 * EN + h * log(numobs)

out <- list(H = H.list, w = w.list, mu = mu.list, psi = psi.list,
            likelihood = likelihood, bic = bic, aic = aic, clc = clc,
            s = s, icl.bic = icl.bic, h = h, ps.y =ps.y)
return(out)
}
