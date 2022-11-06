fix_para <- function(init, init_est) {

  if (any(tolower(init) %in% c('kmeans', 'k-means', 'k')))
    init <- 'kmeans'

  if (any(tolower(init) %in% c('random', 'r')))
    init <- 'random'

  if (any(tolower(init) %in% c('hclass', 'h')))
    init <- 'hclass'

  if (any(tolower(init) %in% c('mclust', 'mclst', 'm')))
    init <- 'mclust'

  if (any(tolower(init_est) == c('factanal', 'factana', 'fact', 'f')))
    init_est <- 'factanal'

  list(init=init, init_est=init_est)
}

