#rm(list=ls())
require(testthat)
y <- mtcars 
layers <- 2 
k <- c(3, 2) 
r <- c(3, 1)
it <- 50 
eps <- 0.001 
seed <- 1 
init <- "kmeans"

model <- deepgmm(y = y, layers = layers, k = k, r = r,
                  it = it, eps = eps, seed = seed, init = init)

expect_that(model, is_a("dgmm"))
#expect_that(model, is_a("emmix"))
expect_named(model, c("H", "w", "mu", "psi", "lik", "bic",
	                    "aic", "clc", "s", "icl.bic", "h", 
                      "k", "r", "numobs"))

# expect_that(g, equals(model$g))
# expect_that(q, equals(model$q))
# expect_that(g, equals(length(model$pivec)))
# expect_that(1, equals(sum(model$pivec)))
# expect_that(p,   equals(nrow(model$A)))
# expect_that(q,   equals(ncol(model$A)))
# expect_that(q,   equals(nrow(model$xi)))
# expect_that(g,   equals(ncol(model$xi)))
# dim_omega <- dim(model$omega)
# expect_that(q, equals(dim_omega[1]))
# expect_that(q, equals(dim_omega[2]))
# expect_that(g, equals(dim_omega[3]))
# expect_that(p, equals(nrow(model$D)))
# expect_that(p, equals(ncol(model$D)))
# expect_that(n, equals(nrow(model$tau)))
# expect_that(g, equals(ncol(model$tau)))
# expect_that(n, equals(nrow(model$Umean)))
# expect_that(q, equals(ncol(model$Umean)))
# expect_that(n, equals(nrow(model$Uclust)))
# expect_that(q, equals(ncol(model$Uclust)))
# dim_U <- dim(model$Uscores)
# expect_that(n, equals(dim_U[1]))
# expect_that(q, equals(dim_U[2]))
# expect_that(g, equals(dim_U[3]))
# expect_that(n, equals(length(model$clust)))

# g <- 1
# q <- 1

# model <- mcfa(Y, g, q, nkmeans = 2, nrandom = 2, tol = 1.e-5)
# expect_that(model, is_a("mcfa"))
# expect_that(model, is_a("emmix"))
# expect_that(g,   equals(model$g))
# expect_that(q,   equals(model$q))
# expect_that(g,  equals(length(model$pivec)))
# expect_that(1,   equals(sum(model$pivec)))
# expect_that(p,   equals(nrow(model$A)))
# expect_that(q,   equals(ncol(model$A)))
# expect_that(q,   equals(nrow(model$xi)))
# expect_that(g,   equals(ncol(model$xi)))
# dim_omega <- dim(model$omega)
# expect_that(q,   equals(dim_omega[1]))
# expect_that(q,   equals(dim_omega[2]))
# expect_that(g,   equals(dim_omega[3]))
# expect_that(p,   equals(nrow(model$D)))
# expect_that(p,   equals(ncol(model$D)))
# expect_that(n,   equals(nrow(model$tau)))
# expect_that(g,   equals(ncol(model$tau)))
# expect_that(n,   equals(nrow(model$Umean)))
# expect_that(q,   equals(ncol(model$Umean)))
# expect_that(n,   equals(nrow(model$Uclust)))
# expect_that(q,   equals(ncol(model$Uclust)))

# dim_U <- dim(model$Uscores)
# expect_that(n,   equals(dim_U[1]))
# expect_that(q,   equals(dim_U[2]))
# expect_that(g,   equals(dim_U[3]))
# expect_that(n,   equals(length(model$clust)))
# context("mctfa")

# q <- 1
# model <- mctfa(Y, g, q, nkmeans = 2, nrandom = 2, tol = 1.e-5)
# expect_named(model, c("g", "q", "pivec", "A", "xi", "omega",
#   "D", "v", "df_update", "logL", "tau", "BIC", "clust", "Uscores", "Uclust", "Umean",
#   "call", "warn_msg"))

# expect_that(model, is_a("mctfa"))
# expect_that(model, is_a("emmix"))
# expect_that(g,   equals(model$g))
# expect_that(q,   equals(model$q))
# expect_that(g,  equals(length(model$pivec)))
# expect_that(1,   equals(sum(model$pivec)))
# expect_that(p,   equals(nrow(model$A)))
# expect_that(q,   equals(ncol(model$A)))
# expect_that(q,   equals(nrow(model$xi)))
# expect_that(g,   equals(ncol(model$xi)))
# dim_omega <- dim(model$omega)
# expect_that(q,   equals(dim_omega[1]))
# expect_that(q,   equals(dim_omega[2]))
# expect_that(g,   equals(dim_omega[3]))
# expect_that(p,   equals(nrow(model$D)))
# expect_that(p,   equals(ncol(model$D)))
# expect_that(n, equals(nrow(model$tau)))
# expect_that(g, equals(ncol(model$tau)))
# expect_that(n, equals(nrow(model$Umean)))
# expect_that(q, equals(ncol(model$Umean)))
# expect_that(n, equals(nrow(model$Uclust)))
# expect_that(q, equals(ncol(model$Uclust)))
# dim_U <- dim(model$Uscores)
# expect_that(n, equals(dim_U[1]))
# expect_that(q, equals(dim_U[2]))
# expect_that(g, equals(dim_U[3]))
# expect_that(n, equals(length(model$clust)))
# expect_that(g, equals(length(model$v)))
