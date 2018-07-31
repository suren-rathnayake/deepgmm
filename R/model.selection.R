## compare different models and return the best one selected according to criterion (BIC or AIC)
## seeds is the different number of seeds to try
## layers can be 1,2 or 3
## k in an integer number of groups

## example: out=model.selection(y,3,2,seeds=1)

model.selection <- function (y, k, layers, seeds = 3, it = 50, eps = 0.001,
                    init = "kmeans", method = "factanal", criterion = "BIC") {

bic.best <- Inf
aic.best <- Inf
p <- dim(y)[2]
pp <- round(p / 2, 0)
ppp <- round(p / 3, 0)
pppp <- round(p / 4, 0)

if (layers == 1) {

  r <- c(1 : pp)

  bic <- array(NA,c(seeds,pp))
  bic.best <- Inf
  aic <- array(NA,c(seeds,pp))
  aic.best <- Inf

  for (i in 1 : seeds) for (rr in 1 : pp) {

      set.seed(i)
      out <- try(deepgmm(y, 1, k, rr, it = it, eps = eps, init = init,
                  method = method))

      if (!is.character(out)) {
        if (criterion=="BIC") if (out$bic<bic.best) {
          out.best <- out
          bic.best <- out$bic
        }
        if (criterion=="AIC")
          if (out$aic<aic.best) {
            out.best <- out
            aic.best <- out$aic
          }

        bic[i, rr] <- out$bic
        aic[i, rr] <- out$aic}
    }

    if (criterion == "BIC")
      index <- which(bic == min(bic, na.rm = TRUE), arr.ind = TRUE)[1, ]

    if (criterion == "AIC")
      index <- which(aic == min(aic, na.rm = TRUE), arr.ind = TRUE)[1, ]

    message("Best Fit: \n")
    cat(paste("Seed=", index[1], " r=", index[2],
              " BIC:", round(out.best$bic, 2), " AIC:",
              round(out.best$aic, 2)))
  }

  if (layers == 2) {

    r <- as.matrix(expand.grid(1 : pp, 1 : ppp))
    r <- r[(r[, 1]) > (r[, 2]), ]
    k <- rbind(c(k, 1), c(k, 2), c(k, 3))

    bic <- array(NA, c(seeds, nrow(k), nrow(r)))
    bic.best <- Inf
    aic <- array(NA, c(seeds, nrow(k), nrow(r)))
    aic.best <- Inf

    for (i in 1 : seeds)
      for (kk in 1 : nrow(k))
        for (rr in 1 : nrow(r)) {
          set.seed(i)
          out <- try(deepgmm(y, 2, k[kk, ], r[rr, ], it = it, eps = eps,
                      init = init, method = method))
          if (!is.character(out)) {
            if (criterion == "BIC")
              if (out$bic < bic.best) {
                out.best <- out
                bic.best <- out$bic
              }
              if (criterion == "AIC")
                if (out$aic < aic.best) {
                  out.best <- out
                  aic.best <- out$aic
                }

              bic[i, kk, rr] <- out$bic
              aic[i, kk, rr] <- out$aic
          }
        }

    if (criterion == "BIC")
      index <- which(bic == min(bic, na.rm = TRUE), arr.ind = TRUE)[1, ]

    if (criterion == "AIC")
      index <- which(aic == min(aic, na.rm = TRUE), arr.ind = TRUE)[1, ]

    message("Best Fit: \n")
    cat(paste("Seed=", index[1], " k=", paste(k[index[2],], collapse=" "),
      " r=", paste(r[index[3], ], collapse = " "),
      " BIC:", round(out.best$bic, 2), " AIC:", round(out.best$aic, 2)))
  }

  if (layers == 3) {

    r <- as.matrix(expand.grid(1 : pp, 1 : ppp, 1 : ppp))
    r <- r[((r[, 1]) > (r[, 2])) & ((r[, 2]) > (r[, 3])), ]

    k <- rbind(c(k, 1, 1), c(k, 2, 1), c(k, 3, 1), c(k, 1, 2),
               c(k, 2, 2), c(k, 3, 2), c(k, 1, 3), c(k, 2, 3),
               c(k, 3, 3))

    bic <- array(NA, c(seeds, nrow(k), nrow(r)))
    bic.best <- Inf
    aic <- array(NA, c(seeds, nrow(k), nrow(r)))
    aic.best <- Inf

    for (i in 1 : seeds)
      for (kk in 1 : nrow(k))
        for (rr in 1 : nrow(r)) {

          set.seed(i)
          out <- try(deepgmm(y, 3, k[kk, ], r[rr, ], it = it, eps = eps,
                      init = init, method = method))
          if (!is.character(out)) {
            if (criterion=="BIC")
              if (out$bic < bic.best) {
                out.best <- out
                bic.best <- out$bic
              }

            if (criterion=="AIC")
            if (out$aic<aic.best) {
              out.best <- out
              aic.best <- out$aic
            }

            bic[i,kk,rr] <- out$bic
            aic[i,kk,rr] <- out$aic
          }
        }

        if (criterion == "BIC")
          index <- which(bic == min(bic, na.rm = TRUE), arr.ind = TRUE)[1, ]

        if (criterion == "AIC")
          index <- which(aic == min(aic, na.rm = TRUE), arr.ind = TRUE)[1, ]

        message("Best Fit: \n")
        cat(paste("Seed=",index[1], " k=", paste(k[index[2], ], collapse=" "),
            " r=", paste(r[index[3], ], collapse=" "),
            " BIC:", round(out.best$bic, 2), " AIC:", round(out.best$aic, 2)))
  }

  out <- list(fit = out.best, bic = bic, aic = aic)
  invisible(out)
}
