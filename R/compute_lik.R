compute.lik <- function(y, numobs, k, mu.list, H.list, psi.list, w.list) {

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

# cc <- 705 - max(pys)
# pys <- pys + cc
# py <- rowSums(exp(pys))
#
# ps.y <- exp(pys) / matrix(py, numobs, tot.k)
# ps.y <- ifelse(is.na(ps.y), 1/k, ps.y)
# py <- exp(-cc) * py

# followning lines were used instead of the lines above
# to avod problmes with near zero pi_i's. 
pys_max <- apply(pys, 1, max)
pys <- sweep(pys, 1, pys_max, '-')
pys <- exp(pys)
ps.y <- sweep(pys, 1, rowSums(pys), '/')
py <- exp(pys_max) * rowSums(pys)

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
