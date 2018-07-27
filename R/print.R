print.dgmm <- function(x, ...) {

  cat("Call:\n")
  print(x$call)

  cat("\n Number of Layers: \n")
  print(x$layers)

  for (lay in 1 : x$layers) {
    cat("\n ---- layer", lay, "-----", "\n")

    cat("\nCoefficients: \n")
    cat("pi_i : ", round(x$w[[lay]], 3), "\n")


	  for (j in 1 : x$k[lay]) {
	    cat("mu_", j, ":\n", sep = "")
	    print(x$mu[[lay]][, j])
	  }

    for (j in 1 : x$k[lay]) {
	    cat("Lambda_", j, ":\n", sep = "")
	    print(x$H[[lay]][j,, ])
	  }

    for (j in 1 : x$k[lay]) {
      cat("diag of Psi_", j, ":\n", sep = "")
      print(diag(x$psi[[lay]][j,, ]))
    }
	}
}

summary.dgmm <- function(object, ...) {

  cat("Call:\n")
  print(object$call)
  summ <- cbind(log_like = object$lik[length(object$lik)],
                BIC = object$bic,
                AIC = object$aic,
                ICL.BIC = object$icl.bic,
                CLC = object$clc
              )
  #cat("\n")
  print(summ)
}
