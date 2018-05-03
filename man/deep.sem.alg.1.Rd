\name{deep.sem.alg.1}
\alias{deep.sem.alg.1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
deep.sem.alg.1(y, numobs, p, r, k, H.list, psi.list, psi.list.inv, mu.list, w.list, z.list, it, eps)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{numobs}{
%%     ~~Describe \code{numobs} here~~
}
  \item{p}{
%%     ~~Describe \code{p} here~~
}
  \item{r}{
%%     ~~Describe \code{r} here~~
}
  \item{k}{
%%     ~~Describe \code{k} here~~
}
  \item{H.list}{
%%     ~~Describe \code{H.list} here~~
}
  \item{psi.list}{
%%     ~~Describe \code{psi.list} here~~
}
  \item{psi.list.inv}{
%%     ~~Describe \code{psi.list.inv} here~~
}
  \item{mu.list}{
%%     ~~Describe \code{mu.list} here~~
}
  \item{w.list}{
%%     ~~Describe \code{w.list} here~~
}
  \item{z.list}{
%%     ~~Describe \code{z.list} here~~
}
  \item{it}{
%%     ~~Describe \code{it} here~~
}
  \item{eps}{
%%     ~~Describe \code{eps} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (y, numobs, p, r, k, H.list, psi.list, psi.list.inv, 
    mu.list, w.list, z.list, it, eps) 
{
    k2 <- 1
    k1 <- k
    r <- r[1]
    lik <- -1e+11
    muf <- matrix(0, p, k1)
    psi <- psi.list[[1]]
    H <- H.list[[1]]
    muf <- mu.list[[1]]
    w1 <- w.list[[1]]
    mu <- array(0, c(1, r))
    sigma <- array(0, c(1, r, r))
    sigma[1, , ] <- diag(r)
    w2 <- 1
    w1 <- matrix(w1)
    w2 <- matrix(w2)
    likelihood <- NULL
    hh <- 0
    ratio <- 1000
    chsi <- array(0, c(k2, k1, r, r))
    roy <- array(0, c(k2, k1, r, numobs))
    py.s1.s2 <- array(0, c(k2, k1, numobs))
    ps1s2.y <- array(0, c(k2, k1, numobs))
    Ezz.y.s1.s2 <- array(0, c(k2, k1, r, r, numobs))
    ps1.y.s2 <- array(0, c(k2, k1, numobs))
    ps2.y.s1 <- array(0, c(k2, k1, numobs))
    sigma.tot <- array(0, c(k2, k1, p, p))
    ps2.y <- array(0, c(k2, numobs))
    ps1.y <- array(0, c(k1, numobs))
    lambda.1 <- rep(1, k1)
    lambda.2 <- rep(1, k2)
    py <- matrix(0, numobs)
    for (i in 1:k2) {
        for (j in 1:k1) {
            sigma.tot[i, j, , ] <- matrix(H[j, , ], ncol = r) \%*\% 
                sigma[i, , ] \%*\% t(matrix(H[j, , ], ncol = r)) + 
                psi[j, , ]
            if (det(as.matrix(sigma.tot[i, j, , ])) < 1e-09) {
                diag(sigma.tot[i, j, , ]) <- diag(sigma.tot[i, 
                  j, , ]) + 0.5
            }
            py.s1.s2[i, j, ] <- dmvnorm(y, muf[, j] + t(matrix(H[j, 
                , ], ncol = r) \%*\% mu[i, ]), as.matrix(sigma.tot[i, 
                j, , ]))
            py.s1.s2 <- ifelse(is.na(py.s1.s2), 1e-10, py.s1.s2)
            py.s1.s2 <- ifelse(py.s1.s2 == 0, 1e-10, py.s1.s2)
            py <- py + w1[j] * w2[i] * py.s1.s2[i, j, ]
        }
    }
    cl <- NULL
    while ((hh < it) & (ratio > eps)) {
        hh <- hh + 1
        Ez.y.s2 <- array(0, c(k2, r, numobs))
        Ez.y.s1 <- array(0, c(k1, r, numobs))
        Ezz.y.s2 <- array(0, c(k2, r, r, numobs))
        Ezz.y.s1 <- array(0, c(k1, r, r, numobs))
        temp <- array(0, c(r, r, numobs))
        for (i in 1:k2) {
            for (j in 1:k1) {
                chsi[i, j, , ] <- ginv(t(H[j, , ]) \%*\% ginv(psi[j, 
                  , ]) \%*\% H[j, , ] + ginv(sigma[i, , ]))
                roy[i, j, , ] <- chsi[i, j, , ] \%*\% (t(H[j, , 
                  ]) \%*\% ginv(psi[j, , ]) \%*\% t(y - t(matrix(muf[, 
                  j], p, numobs))) + matrix(ginv(sigma[i, , ]) \%*\% 
                  mu[i, ], r, numobs))
                roy[i, j, , ] <- ifelse(is.na(roy[i, j, , ]), 
                  0, roy[i, j, , ])
                if (r > 1) {
                  for (h in 1:numobs) {
                    temp[, , h] <- (roy[i, j, , h]) \%*\% t(roy[i, 
                      j, , h])
                  }
                  temp2 <- array(chsi[i, j, , ], c(r, r, numobs))
                  Ezz.y.s1.s2[i, j, , , ] <- temp + temp2
                }
                else {
                  Ezz.y.s1.s2[i, j, r, r, ] <- roy[i, j, , ]^2 + 
                    rep(chsi[i, j, , ], numobs)
                }
            }
        }
        for (i in 1:k2) {
            ps2.y[i, ] <- (w2[i] * (t(w1) \%*\% py.s1.s2[i, , ]))/t(py)
            ps2.y[i, ] <- ifelse(is.na(ps2.y[i, ]), mean(ps2.y[i, 
                ], na.rm = TRUE), ps2.y[i, ])
            ps2.y.s1[i, , ] <- w2[i] * py.s1.s2[i, , ]
        }
        A <- apply(ps2.y.s1, c(2, 3), sum)
        A <- array(A, c(k1, numobs, k2))
        A <- aperm(A, c(3, 1, 2))
        A <- ifelse(A == 0, 1e-06, A)
        ps2.y.s1 <- ps2.y.s1/A
        for (j in 1:k1) {
            ps1.y[j, ] <- (w1[j] * (t(w2) \%*\% py.s1.s2[, j, ]))/t(py)
            ps1.y[j, ] <- ifelse(is.na(ps1.y[j, ]), mean(ps1.y[j, 
                ], na.rm = TRUE), ps1.y[j, ])
            ps1.y.s2[, j, ] <- w1[j] * py.s1.s2[, j, ]
        }
        A <- apply(ps1.y.s2, c(1, 3), sum)
        A <- array(A, c(k2, numobs, k1))
        A <- aperm(A, c(1, 3, 2))
        A <- ifelse(A == 0, 1e-06, A)
        ps1.y.s2 <- ps1.y.s2/A
        for (h1 in 1:r) {
            for (j in 1:k1) {
                Ez.y.s2[, h1, ] <- Ez.y.s2[, h1, ] + roy[, j, 
                  h1, ] * ps1.y.s2[, j, ]
                for (h2 in 1:r) {
                  Ezz.y.s2[, h1, h2, ] <- Ezz.y.s2[, h1, h2, 
                    ] + Ezz.y.s1.s2[, j, h1, h2, ] * ps1.y.s2[, 
                    j, ]
                }
            }
            for (i in 1:k2) {
                Ez.y.s1[, h1, ] <- Ez.y.s1[, h1, ] + roy[i, , 
                  h1, ] * ps2.y.s1[i, , ]
                for (h2 in 1:r) {
                  Ezz.y.s1[, h1, h2, ] <- Ezz.y.s1[, h1, h2, 
                    ] + Ezz.y.s1.s2[i, , h1, h2, ] * ps2.y.s1[i, 
                    , ]
                }
            }
        }
        temp <- aperm(array(ps1.y, c(k1, numobs, r)), c(1, 3, 
            2))
        Ez.y <- Ez.y.s1 * temp
        Ez.y <- t(apply(Ez.y, c(2, 3), sum))
        w2 <- rowMeans(ps2.y)
        w1 <- rowMeans(ps1.y)
        for (i in 1:k2) {
            Ez.y.s2[i, , ] <- ifelse(is.na(Ez.y.s2[i, , ]), rowMeans(matrix(Ez.y.s2[i, 
                , ], ncol = numobs), na.rm = TRUE), Ez.y.s2[i, 
                , ])
            if (r > 1) {
                Ezz.y.s2[i, , , ] <- ifelse(is.na(Ezz.y.s2[i, 
                  , , ]), apply(Ezz.y.s2[i, , , ], c(1, 2), mean, 
                  na.rm = TRUE), Ezz.y.s2[i, , , ])
            }
            if (r == 1) {
                Ezz.y.s2[i, , , ] <- ifelse(is.na(Ezz.y.s2[i, 
                  , , ]), mean(Ezz.y.s2[i, , , ], na.rm = TRUE), 
                  Ezz.y.s2[i, , , ])
            }
            mu[i, ] <- t((Ez.y.s2[i, , ] \%*\% ps2.y[i, ])/sum(ps2.y[i, 
                ]))
            sigma[i, , ] <- apply(((Ezz.y.s2[i, , , ] - array(mu[i, 
                ] \%*\% t(mu[i, ]), c(r, r, numobs))) * (aperm(array(matrix(ps2.y[i, 
                ]), c(numobs, r, r)), c(2, 3, 1)))), 1, rowSums)/sum(ps2.y[i, 
                ])
        }
        for (j in 1:k1) {
            Ez.y.s1[j, , ] <- ifelse(is.na(Ez.y.s1[j, , ]), rowMeans(matrix(Ez.y.s1[j, 
                , ], ncol = numobs), na.rm = TRUE), Ez.y.s1[j, 
                , ])
            if (r > 1) {
                Ezz.y.s1[j, , , ] <- ifelse(is.na(Ezz.y.s1[j, 
                  , , ]), apply(Ezz.y.s1[j, , , ], c(1, 2), mean, 
                  na.rm = TRUE), Ezz.y.s1[j, , , ])
            }
            if (r == 1) {
                Ezz.y.s1[j, , , ] <- ifelse(is.na(Ezz.y.s1[j, 
                  , , ]), mean(Ezz.y.s1[j, , , ], na.rm = TRUE), 
                  Ezz.y.s1[j, , , ])
            }
            if (r > 1) {
                EEzz.y.s1 <- apply((Ezz.y.s1[j, , , ] * (aperm(array(t(t(ps1.y[j, 
                  ])), c(numobs, r, r)), c(2, 3, 1)))), 1, rowSums)/sum(ps1.y[j, 
                  ])
            }
            if (r == 1) {
                EEzz.y.s1 <- Ezz.y.s1[j, , , ] \%*\% (ps1.y[j, 
                  ])/sum(ps1.y[j, ])
            }
            H[j, , ] <- (t((y - t(matrix(muf[, j], p, numobs))) * 
                matrix(ps1.y[j, ], numobs, p)) \%*\% (t(matrix(Ez.y.s1[j, 
                , ], ncol = numobs)) * matrix(ps1.y[j, ], numobs, 
                r))) \%*\% ginv(EEzz.y.s1)/sum(ps1.y[j, ])
            H[j, , ] <- ifelse(is.na(H[j, , ]), 0.1, H[j, , ])
            psi[j, , ] <- (t((y - t(matrix(muf[, j], p, numobs))) * 
                matrix(ps1.y[j, ], numobs, p)) \%*\% (matrix((y - 
                t(matrix(muf[, j], p, numobs))) * matrix(ps1.y[j, 
                ], numobs, p), ncol = p)) - t((y - t(matrix(muf[, 
                j], p, numobs))) * matrix(ps1.y[j, ], numobs, 
                p)) \%*\% (t(matrix(Ez.y.s1[j, , ], ncol = numobs)) * 
                matrix(ps1.y[j, ], numobs, r)) \%*\% t(H[j, , ]))
            muf[, j] <- colSums(matrix(ps1.y[j, ], numobs, p) * 
                (y - t(matrix(H[j, , ], ncol = r) \%*\% Ez.y.s1[j, 
                  , ])))/sum(ps1.y[j, ])
        }
        for (j in 1:k1) {
            psi[j, , ] <- psi[j, , , drop = FALSE]/sum(ps1.y[j, 
                ])
            if (p > 1) {
                psi[j, , ] <- diag(diag(psi[j, , ]))
                psi[j, , ] <- ifelse(is.na(psi[j, , ]), 1, psi[j, 
                  , ])
            }
        }
        py <- matrix(0, numobs)
        for (i in 1:k2) {
            for (j in 1:k1) {
                sigma.tot[i, j, , ] <- matrix(H[j, , ], ncol = r) \%*\% 
                  sigma[i, , ] \%*\% t(matrix(H[j, , ], ncol = r)) + 
                  psi[j, , ]
                if (det(as.matrix(sigma.tot[i, j, , ])) < 1e-09) {
                  diag(sigma.tot[i, j, , ]) <- diag(sigma.tot[i, 
                    j, , ]) + 0.5
                }
                py.s1.s2[i, j, ] <- dmvnorm(y, muf[, j] + t(matrix(H[j, 
                  , ], ncol = r) \%*\% mu[i, ]), as.matrix(sigma.tot[i, 
                  j, , ]))
                py.s1.s2 <- ifelse(is.na(py.s1.s2), 1e-07, py.s1.s2)
                ps1s2.y[i, j, ] <- w1[j] * w2[i] * py.s1.s2[i, 
                  j, ]
                py <- py + ps1s2.y[i, j, ]
            }
        }
        ps1s2.y <- ps1s2.y/aperm(array(py, c(numobs, k2, k1)), 
            c(2, 3, 1))
        temp <- sum(log(py))
        likelihood <- c(likelihood, temp)
        ratio <- (temp - lik)/abs(lik)
        if (hh < 5) {
            ratio <- 2 * eps
        }
        lik <- temp
    }
    if (k1 > 1) {
        s1 <- apply(ps1.y, 2, order)[k1, ]
    }
    else {
        s1 <- rep(k1, numobs)
    }
    likelihood <- matrix(likelihood[!likelihood == 0])
    h1 <- (k1 - 1) + (p * r) * k1 + p * k1 - (k1 * r * (r - 1)/2) + 
        (k1 - 1) * p
    h2 <- (k2 - 1) + r * (k2 - 1) + (r * (r + 1)/2) * (k2 - 1)
    h <- h1 + h2
    lik <- likelihood[length(likelihood)]
    bic <- -2 * lik + h * log(numobs)
    aic <- -2 * lik + 2 * h
    EN <- entr(matrix(ps1s2.y, k2 * k1, numobs))
    clc <- -2 * lik + 2 * EN
    icl.bic <- -2 * lik + 2 * EN + h * log(numobs)
    out <- list(H = H, w = w1, mu = muf, psi = psi, likelihood = likelihood, 
        bic = bic, aic = aic, clc = clc, s = s1, icl.bic = icl.bic, 
        h = h, ps.y = ps.y)
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
