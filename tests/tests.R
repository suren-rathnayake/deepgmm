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

