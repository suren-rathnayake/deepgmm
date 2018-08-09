deepgmm <- function(y, layers, k, r = rep(1, layers),
            it = 50, eps = 0.001, init = 'kmeans', method = "factanal") {

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

    if (init == 'kmeans') {
      if (k[i] > 1) {
        s <- kmeans(data, k[i], iter.max = 100, nstart = 30,
               algorithm = "Hartigan-Wong")$cluster
      } else {
        s <- rep(1, numobs)
      }
    }

    if (init == 'hclass') {
      if (k[i] > 1) {
        s <- cutree(hclust(dist(y), "ward.D2"), k[i])
      } else {
        s <- rep(1, numobs)
      }
    }

    if (init == 'random') {
      if (k[i] > 1) {
        s <- sample(1 : k[i], numobs, replace = TRUE)
      } else {
        s <- rep(1, numobs)
      }
    }

    for  (j in 1 : k[i]) {
      if ((table(s)[j]) < 2) {
        s[sample(1 : numobs, 2, replace = FALSE)] <- j
      }
    }

    psi <- psi.inv <- array(0, c(k[i], r[i], r[i]))
    H <- array(0, c(k[i], r[i], r[i + 1]))
    mu <- matrix(0, r[i], k[i])
	  if (method == "factanal") {

      z <- matrix(NA, nrow = numobs, ncol = r[i + 1])
	    for (j in 1 : k[i]) {

	    	indices <- which(s == j)
	      stima <- try(factanal(data[indices, ], r[i + 1], rotation = "none",
	                  scores = "Bartlett"), silent = TRUE)

	      if (is.character(stima)) {

          psi[j,, ] <- 0.1 * diag(r[i])
	        psi.inv[j,, ] <- diag(r[i])
	        H[j,, ] <- matrix(runif(r[i] * r[i + 1]), r[i], r[i + 1])

	        zt <- try(princomp(data[indices, ])$scores[, 1 : r[i + 1]],
                    silent = TRUE)
	        if (!is.character(zt)) {
	          zt <- matrix(zt, ncol = r[i+1])
	        } else {
	          zt <- matrix(data[indices, sample(1 : r[i + 1])], ncol = r[i + 1])
	        }
	        z[indices, ] <- zt
	      } else {

	        psi[j,, ] <- diag(stima$uniq)
	        H[j,, ] <- stima$load
	        psi.inv[j,, ] <- diag(1/stima$uniq)
	        z[indices, ] <- stima$scores
	      }

	      mu[, j] <- colMeans(data[indices,, drop = FALSE])
	    }

    } else {

      z <- matrix(NA, nrow = numobs, ncol = r[i + 1])
			for (j in 1 : k[i]) {

				q <- r[i + 1]
			  indices <- which(s == j)
			  mu[, j] <- colMeans(data[indices,, drop = FALSE])
			  Si <- cov(data[indices, ])
			  psi[j,, ] <-  diag(diag(Si))
			  Di.sqrt <- diag(sqrt(diag(diag(diag(Si)))))
			  inv.Di.sqrt <- diag(1 / diag(Di.sqrt))
			  eig.list <- eigen(inv.Di.sqrt %*% Si %*% inv.Di.sqrt)
			  eigH <- eig.list$vectors
			  sort.lambda <- sort(eig.list$values, decreasing = TRUE,
			                                       index.return = TRUE)
			  lambda <- sort.lambda$x
			  ix.lambda   <- sort.lambda$ix
			  sigma2 <- mean(lambda[(q + 1) : ncol(data)])
			  if (q == 1) {
			    H[j,, ] <- Di.sqrt %*% eigH[, ix.lambda[1 : q]] %*%
			                    diag((lambda[1 : q] - sigma2), q)
			    z[indices, ] <-  sweep(data[indices,, drop = FALSE], 2,
			   	                   mu[, j, drop = FALSE], '-') %*%
			                        t(1 / (t(H[j,, ]) %*% H[j,, ] +
			                        	diag(sigma2, q)) %*% t(H[j,, ]))
			  } else {
			    H[j,, ] <- Di.sqrt %*% eigH[, ix.lambda[1 : q]] %*%
			                    diag((lambda[1 : q] - sigma2))

	    		z[indices, ] <-  sweep(data[indices,, drop = FALSE], 2,
	    			                   mu[, j, drop = FALSE], '-') %*%
			                         t(chol.inv(t(H[j,, ]) %*% H[j,, ] +
			                         	diag(sigma2, q)) %*% t(H[j,, ]))
			  }
			}
		}

    w <- matrix(table(s) / numobs)
    lst$w[i] <- list(w)
    lst$H[i] <- list(H)
    lst$mu[i] <- list(mu)
    lst$psi[i] <- list(psi)
    lst$psi.inv[i] <- list(psi.inv)
  }

  ##############################################################################
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

  # if (! class(s) %in% "matrix")
  #   s <- matrix(s, nrow = numobs)

  out$lik <- out$likelihood
  output <- out[c("H", "w", "mu", "psi", "lik", "bic", "aic", "clc",
                  "icl.bic", "s", "h")]
  output <- c(output, list(k = k, r = r[-1], numobs = numobs, layers = layers))
  output$call <- match.call()
  class(output) <- "dgmm"

  message("Estimation Details:")
  cat(paste("Log-Likelihood", round(output$lik[length(output$lik)], 2), "BIC:",
        round(output$bic, 2), "AIC:", round(output$aic, 2), "\n"))
  cat("\n")
  invisible(output)
}
