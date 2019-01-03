deepgmm <- function(y, layers, k, r,
            it = 250, eps = 0.001, init = 'kmeans', method = "factanal") {

  if (any(tolower(init) == c('kmeans', 'k-means', 'k')))
    init <- 'kmeans'
  if (any(tolower(init) == c('random', 'r')))
    init <- 'random'
  if (any(tolower(init) == c('hclass', 'h')))
    init <- 'hclass'
  if (any(tolower(method) == c('factanal', 'factana', 'fact', 'f')))
    method <- "factanal"
  if (class(y) == "data.frame")
  	y <- as.matrix(y)

  # check arguments
  tmp <- valid_args(Y = y, layers = layers, k = k, r = r, it = it,
                    eps = eps, init = init)
  numobs <- nrow(y)
  p <- ncol(y)
  r <- c(p, r)

  # Initialing parameters
  lst <- list(w = list(), H = list(), mu = list(), psi = list(),
                                                   psi.inv = list())
  for (i in 1 : layers) {

    if (i == 1) {
      data <- y
    } else {
      data <- z[, 1 : r[i], drop = FALSE]
    }

    # provide initial parititioning of the observations
    s <- initial_clustering(data, k, i, init)

    # in case if one of the groups is small
    for  (j in 1 : k[i]) {
      if ((table(s)[j]) < 2) {
        s[sample(1 : numobs, 2, replace = FALSE)] <- j
      }
    }

    psi <- psi.inv <- array(0, c(k[i], r[i], r[i]))
    H <- array(0, c(k[i], r[i], r[i + 1]))
    mu <- matrix(0, r[i], k[i])

	  if (method == "factanal") {
        i_lst <- factanal_para(data, s, k, r, i, numobs)

      lst$w[i] <- list(i_lst$w)
      lst$H[i] <- list(i_lst$H)
      lst$mu[i] <- list(i_lst$mu)
      lst$psi[i] <- list(i_lst$psi)
      lst$psi.inv[i] <- list(i_lst$psi.inv)
      z <- i_lst$z

    } else {
    # initialize parameters using probabilistic principal component analysis

      if (method != "ppca")
        stop("method has to be either `factanal` or `ppca`")

      # z <- matrix(NA, nrow = numobs, ncol = r[i + 1])
			# for (j in 1 : k[i]) {
      #
			# 	q <- r[i + 1]
			#   indices <- which(s == j)
			#   mu[, j] <- colMeans(data[indices,, drop = FALSE])
			#   Si <- cov(data[indices, ])
			#   psi[j,, ] <-  diag(diag(Si))
			#   Di.sqrt <- diag(sqrt(diag(diag(diag(Si)))))
			#   inv.Di.sqrt <- diag(1 / diag(Di.sqrt))
			#   eig.list <- eigen(inv.Di.sqrt %*% Si %*% inv.Di.sqrt)
			#   eigH <- eig.list$vectors
			#   sort.lambda <- sort(eig.list$values, decreasing = TRUE,
			#                                        index.return = TRUE)
			#   lambda <- sort.lambda$x
			#   ix.lambda   <- sort.lambda$ix
			#   sigma2 <- mean(lambda[(q + 1) : ncol(data)])
			#   if (q == 1) {
			#     H[j,, ] <- Di.sqrt %*% eigH[, ix.lambda[1 : q]] %*%
			#                     diag((lambda[1 : q] - sigma2), q)
			#     z[indices, ] <-  sweep(data[indices,, drop = FALSE], 2,
			#    	                   mu[, j, drop = FALSE], '-') %*%
			#                         t(1 / (t(H[j,, ]) %*% H[j,, ] +
			#                         	diag(sigma2, q)) %*% t(H[j,, ]))
			#   } else {
			#     H[j,, ] <- Di.sqrt %*% eigH[, ix.lambda[1 : q]] %*%
			#                     diag((lambda[1 : q] - sigma2))
      #
	    # 		z[indices, ] <-  sweep(data[indices,, drop = FALSE], 2,
	    # 			                   mu[, j, drop = FALSE], '-') %*%
			#                          t(chol.inv(t(H[j,, ]) %*% H[j,, ] +
			#                          	diag(sigma2, q)) %*% t(H[j,, ]))
			#   }
			# }
      #
      # w <- matrix(table(s) / numobs)
      # lst$w[i] <- list(w)
      # lst$H[i] <- list(H)
      # lst$mu[i] <- list(mu)
      # lst$psi[i] <- list(psi)
      # lst$psi.inv[i] <- list(psi.inv)
      #

      i_lst <- ppca_para(data, s, k, r, i, numobs)
      lst$w[i] <- list(i_lst$w)
      lst$H[i] <- list(i_lst$H)
      lst$mu[i] <- list(i_lst$mu)
      lst$psi[i] <- list(i_lst$psi)
      lst$psi.inv[i] <- list(i_lst$psi.inv)
      z <- i_lst$z
    }


  }

  ############################################################################
  if (layers == 1) {
    out <- deep.sem.alg.1(y, numobs, p, r[2], k, lst$H, lst$psi,
                        lst$psi.inv, lst$mu, lst$w, it, eps)

  }
  if (layers == 2) {
    out <- deep.sem.alg.2(y, numobs, p, r, k, lst$H, lst$psi,
                        lst$psi.inv, lst$mu, lst$w, it, eps)
  }
  if (layers == 3) {
    out <- deep.sem.alg.3(y, numobs, p, r, k, lst$H, lst$psi,
                        lst$psi.inv, lst$mu, lst$w, it, eps)
  }

  out$lik <- out$likelihood
  output <- out[c("H", "w", "mu", "psi", "lik", "bic", "aic", "clc",
                  "icl_bic", "s", "h")]
  output <- c(output, list(k = k, r = r[-1], numobs = numobs, layers = layers))
  output$call <- match.call()
  class(output) <- "dgmm"

  # message("Estimation Details:")
  # cat(paste("Log-Likelihood", round(output$lik[length(output$lik)], 2), "BIC:",
  #       round(output$bic, 2), "AIC:", round(output$aic, 2), "\n"))
  # cat("\n")
  invisible(output)
}
