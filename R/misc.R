chol.inv <- function(x, ...) {
  C <- chol(x)
  inv_x <- chol2inv(C)
  return(inv_x)
}

makeSymm <- function(m) {
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  return(m)

  # pmean <- function(x, y) (x + y) / 2
  # m[] <- pmean(m, matrix(m, nrow(m), byrow=TRUE))
  # m
}

ma <- function(x, n = 5) {
  filter(x, rep(1/n, n), sides = 1)
}

entr <- function(z) {
  numobs <- nrow(z)
  numg <- ncol(z)
  temp <- 0
  z <- ifelse(z == 0, z + .Machine$double.eps, z)
  for (i in 1 : numg)  {
    for (j in 1 : numobs) {
      temp <- temp + (z[j, i] * log(z[j, i]))
    }
  }
  return(-temp)
}

