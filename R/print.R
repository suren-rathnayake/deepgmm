print.dgmm <- function(x, ...) {

  cat("Call:\n")
  print(x$call)

  cat("\n Number of Layers: \n")
  print(x$layers)

  for (lay in 1 : x$layers) {
    cat(" ---- layer", lay, "-----", "\n")

    cat("\nCoefficients: \n")
    cat("pi_i : ", round(x$w[[lay]], 3), "\n")


	  for(j in 1 : x$k[lay]) {
	    cat("mu_", j, ":\n")
	    print(x$mu[[lay]][, j])
	  }


	}

}

summary.dgmm <- function(object, ...) {

  cat("Call:\n")
  print(object$call)
  summ <- cbind(log_like = object$logL,
                BIC = object$bic,
                AIC = object$aic,
                ICL.BIC = object$cls,
                CIC = object$icl.bic
              )
  cat("\n")
  print(summ)
}
