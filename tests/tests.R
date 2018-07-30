#rm(list=ls())
require(testthat)
library(deepgmm)
y <- mtcars
layers <- 2
k <- c(3, 2)
r <- c(3, 1)
it <- 50
eps <- 0.001
seed <- 1
init <- "random"

set.seed(seed)
y <- scale(y)
model <- deepgmm(y = y, layers = layers, k = k, r = r,
                  it = it, eps = eps, init = init)

expect_that(model, is_a("dgmm"))
#expect_that(model, is_a("emmix"))
expect_named(model, c("H", "w", "mu", "psi", "lik", "bic",
	                    "aic", "clc", "s", "icl.bic", "h",
                      "k", "r", "numobs", "layers",  "call"))

n <- nrow(y)
expect_that(layers, equals(model$layers))
expect_that(n, equals(model$numobs))
expect_that(k, equals(model$k))
expect_that(r, equals(model$r))

for (j in 1 : layers) {
  expect_that(n, equals(length(model$s[, j])))
  expect_that(k[j], equals(length(model$w[[j]])))

}
# expect_that(g, equals(length(model$pivec)))
# expect_that(1, equals(sum(model$pivec)))
#
# expect_that(p, equals(nrow(model$mu)))
# expect_that(g, equals(ncol(model$mu)))
#
# dim_D <- dim(model$D)
# expect_that(p, equals(dim_D[1]))
# expect_that(p, equals(dim_D[2]))
#
# dim_sigma <- dim(model$B)
# expect_that(p, equals(dim_sigma[1]))
# expect_that(q, equals(dim_sigma[2]))
# expect_that(g, equals(dim_sigma[3]))
#
# expect_that(n, equals(nrow(model$tau)))
# expect_that(g, equals(ncol(model$tau)))
#
# expect_that(n, equals(nrow(model$Umean)))
# expect_that(q, equals(ncol(model$Umean)))
#
# expect_that(n, equals(nrow(model$Uclust)))
# expect_that(q, equals(ncol(model$Uclust)))
