deepgmm <- function(y, layers, k, r,
            it = 250, eps = 0.001, init = 'kmeans', init_est = 'factanal',
            seed = NULL, scale = TRUE) {

  if (any(class(y) %in% 'data.frame'))
  	y <- as.matrix(y)

  if (!is.null(seed)) {

    if(!is.numeric(seed)) {
      stop("The value of seed must be an integer")
    }

    set.seed(seed)
  }

  if (scale)
    y <- scale(y)

  if (any(tolower(init) %in% c('kmeans', 'k-means', 'k')))
    init <- 'kmeans'
  if (any(tolower(init) %in% c('random', 'r')))
    init <- 'random'
  if (any(tolower(init) %in% c('hclass', 'h')))
    init <- 'hclass'
  if (any(tolower(init_est) == c('factanal', 'factana', 'fact', 'f')))
    init_est <- 'factanal'


  # check arguments
  tmp <- valid_args(Y = y, layers = layers, k = k, r = r, it = it,
                    eps = eps, init = init)
  numobs <- nrow(y)
  p <- ncol(y)
  r <- c(p, r)

  # Initialing parameters
  lst <- list(w = list(), H = list(), mu = list(), psi = list(),
                                                   psi.inv = list())
for (i in 1:layers) {

    if (i == 1) {

      data <- y
    } else {

      data <- z[, 1 : r[i], drop = FALSE]
    }

    # provide initial parititioning of the observations
    s <- initial_clustering(data, k, i, init)

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

	  if (init_est == "factanal") {
    # initialize parameters using factor analysis of covariance matrix

      i_lst <- factanal_para(data, s, k, r, i, numobs)

      lst$w[i] <- list(i_lst$w)
      lst$H[i] <- list(i_lst$H)
      lst$mu[i] <- list(i_lst$mu) 
      lst$psi[i] <- list(i_lst$psi)
      lst$psi.inv[i] <- list(i_lst$psi.inv)
      z <- i_lst$z

    } else {
    # initialize parameters using probabilistic principal component analysis

      if (init_est != "ppca")
        stop("init_est has to be either `factanal` or `ppca`")

      i_lst <- ppca_para(data, s, k, r, i, numobs)

      lst$w[i] <- list(i_lst$w)
      lst$H[i] <- list(i_lst$H)
      lst$mu[i] <- list(i_lst$mu)
      lst$psi[i] <- list(i_lst$psi)
      lst$psi.inv[i] <- list(i_lst$psi.inv)
      z <- i_lst$z
    }
  }

  ##########################################################
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
  output <- c(output, list(k = k, r = r[-1], numobs = numobs, layers = layers,
                           seed = seed))
  output$call <- match.call()
  class(output) <- "dgmm"

  invisible(output)
}
