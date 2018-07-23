deepgmm <- function(y, layers, k, r = rep(1, layers),
                  it = 50, eps = 0.001, init = 'kmeans', method = "factanal") {

  if (any(tolower(init) == c('kmeans', 'k-means', 'k')))
    init <- 'kmeans'

  if (any(tolower(init) == c('random', 'r')))
    init <- 'random'

  if (any(tolower(init) == c('hclass', 'h')))
    init <- 'hclass'

  if (class(y) == "data.frame") 
  	y <- as.matrix(y)

  # check arguments
  tmp <- valid_args(Y = y, layers = layers, k = k, r = r, it = it, eps = eps, 
						 init = init)

  # ptm <- proc.time()
  # a <- paste("Fit the model: seed=", seed, " k = (",paste(k, collapse=" "),
  # 	   ") r = (",paste(r, collapse = " "), ")\n", sep = "")
  # cat(a)
  # set.seed(seed)
  # y <- scale(y)

  numobs <- nrow(y)
  p <- ncol(y)
  r <- c(p, r) #### *

  #init
  w.list <- NULL
  H.list <- NULL
  psi.list <- NULL
  psi.list.inv <- NULL
  mu.list <- NULL
  z.list <- NULL

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
        #s <- hclass(hc(modelName = "VVV", data = y), k[i])
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
    H <- array(0, c(k[i], r[i], r[i+1]))
    mu <- matrix(0, r[i], k[i])

     z <- NULL
    #z <- matrix(NA, nrow = numobs, ncol = r[i + 1])

	  if (method == "factanal") {

	    for (j in 1 : k[i]) {

	    	indices <- which(s == j)

	      stima <- try(factanal(data[s == j, ], r[i + 1], rotation = "none",
	                  scores = "Bartlett"), silent = TRUE)

	      if (is.character(stima)) {

	        psi[j,, ] <- 0.1 * diag(r[i])
	        psi.inv[j,, ] <- diag(r[i])
	        H[j,,] <- matrix(runif(r[i] * r[i+1]), r[i], r[i+1])
	        zt <- try(princomp(data[s == j, ])$scores[, 1 : r[i+1]], silent = TRUE)

	        if (!is.character(zt)) {

	          zt <- matrix(zt, ncol = r[i+1])
	        }

	        if (is.character(zt)) {

	          zt <- matrix(data[s == j, sample(1 : r[i+1])], ncol = r[i+1])
	        }
	        z <- rbind(z, zt)
	        #z[indices, ] <- zt
	      }

	      if (!is.character(stima)) {

	        psi[j,, ] <- diag(stima$uniq)
	        H[j,, ] <- stima$load
	        psi.inv[j,,] <- diag(1/stima$uniq)
	        z <- rbind(z, stima$scores)
	        #z[indices, ] <- stima$scores
	      }
	      
	      mu[, j] <- colMeans(data[s == j,, drop = FALSE])
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
			  # eig.list <- try(eigen(inv.Di.sqrt %*% Si %*% inv.Di.sqrt), TRUE)
			  # if (class(eig.list) == "try-error")
			  #   break
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
			                        	diag(lambda[1 : q], q)) %*% t(H[j,, ]))                    
			  } else {
			    H[j,, ] <- Di.sqrt %*% eigH[, ix.lambda[1 : q]] %*%
			                    diag((lambda[1 : q] - sigma2))

	    		z[indices, ] <-  sweep(data[indices,, drop = FALSE], 2, 
	    			                   mu[, j, drop = FALSE], '-') %*%
			                         t(chol.inv(t(H[j,, ]) %*% H[j,, ] + 
			                         	diag(mean(lambda[1 : q]), q) %*% t(H[j,, ]))                    
			  }
			}
		}

    w <- matrix(table(s) / numobs)
    w.list[i] <- list(w)
    H.list[i] <- list(H)
    mu.list[i] <- list(mu)
    psi.list[i] <- list(psi)
    psi.list.inv[i] <- list(psi.inv)
    z.list[i] <- list(aperm(array(z[, 1 : r[i+1]], c(numobs, r[i+1], k[i])),
                     c(3, 1, 2)))
  }

  ##############################################################################
  if (layers == 1) {
    out <- deep.sem.alg.1(y, numobs, p, r[2], k, H.list, psi.list,
                        psi.list.inv, mu.list, w.list, z.list, it, eps) 
  }
  if (layers == 2) {
    out <- deep.sem.alg.2(y, numobs, p, r, k, H.list, psi.list,
                        psi.list.inv, mu.list, w.list, z.list, it, eps)
  }
  if (layers == 3) {
    out <- deep.sem.alg.3(y, numobs, p, r, k, H.list, psi.list,
                        psi.list.inv, mu.list, w.list, z.list, it, eps)
  }

  H <- out$H
  w <- out$w
  mu <- out$mu
  psi <- out$psi
  lik <- out$likelihood
  s <- out$s
  bic <- out$bic
  h <- out$h
  aic <- out$aic
  icl.bic <- out$icl.bic
  clc <- out$clc

  output <- list (H = H, w = w, mu = mu, psi = psi, lik = lik,
                 bic = bic, aic = aic, clc = clc, s = s, icl.bic = icl.bic,
                 h = h, k = k, r = r, numobs = numobs, layers = layers)

  output$call <- match.call()
  class(output) <- "dgmm"

  message("Estimation Details:")
  cat(paste("Log-Likelihood", round(lik[length(lik)], 2), "BIC:",
        round(bic, 2), "AIC:", round(aic, 2), "\n"))
  cat("\n")
  invisible(output)
}