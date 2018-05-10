print.dgmm <- function(x, ...) {

  cat("Call:\n")
  print(x$call)

  cat("\n Number of Layers: \n")
  print(x$layers)

  for (lay in 1 : layers) {
    cat(" ---- layer", lay, "-----", "\n")

    cat("\nCoefficients: \n")
    cat("pi_i : ", round(x$w[[lay]], 3), "\n")

  
	  for(j in 1 : x$k[lay]) {
	    cat("mu_", j, ":\n")
	    print(x$mu[[lay]][, j])
	  }

	  
	}

}
