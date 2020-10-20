initial_clustering <- function(data, k, i, init = 'random') {

  numobs <- nrow(data)

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
      s <- cutree(hclust(dist(data), "ward.D2"), k[i])
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

  if (init == 'mclust') {
   if (k[i] > 1) {
      s <- Mclust(data, k[i], verbose=FALSE)$classification
    } else { 
      s <- rep(1, numobs)
    }
  }

  s
}
