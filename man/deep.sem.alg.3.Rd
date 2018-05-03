\name{deep.sem.alg.3}
\alias{deep.sem.alg.3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
deep.sem.alg.3(y, numobs, p, r, k, H.list, psi.list, psi.list.inv, mu.list, w.list, z.list, it, eps)
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
    likelihood <- NULL
    hh <- 0
    ratio <- 1000
    layers = length(k)
    out = compute.lik(y, numobs, k, mu.list, H.list, psi.list, 
        w.list)
    py = out$py
    ps.y = out$ps.y
    ps.y.list = out$ps.y.list
    k.comb = out$k.comb
    s = out$s
    tot.k = prod(k)
    temp = sum(log(py))
    likelihood <- c(likelihood, temp)
    while ((hh < it) & (ratio > eps)) {
        hh <- hh + 1
        l = 1
        yy = y
        z2 <- z2.one <- array(0, c(numobs, r[l + 1], k[l], k[l + 
            1], k[l + 2]))
        zz2 <- array(0, c(numobs, r[l + 1], r[l + 1], k[l], k[l + 
            1], k[l + 2]))
        z <- z.one <- array(0, c(numobs, r[l + 1], k[l]))
        zz <- array(0, c(numobs, r[l + 1], r[l + 1], k[l]))
        for (p1 in 1:k[l]) for (p2 in 1:k[l + 1]) for (p3 in 1:k[l + 
            2]) {
            sigma.tilde.inv = ginv(H.list[[l + 1]][p2, , ] \%*\% 
                (H.list[[l + 2]][p3, , ] \%*\% t(H.list[[l + 2]][p3, 
                  , ]) + psi.list[[l + 2]][p3, , ]) \%*\% t(H.list[[l + 
                1]][p2, , ]) + psi.list[[l + 1]][p2, , ])
            A = sigma.tilde.inv + t(H.list[[l]][p1, , ]) \%*\% 
                (psi.list.inv[[l]][p1, , ]) \%*\% (H.list[[l]][p1, 
                , ])
            mu.tilde = matrix(mu.list[[l + 1]][, p2] + H.list[[l + 
                1]][p2, , ] \%*\% (mu.list[[l + 2]][, p3]), r[l + 
                1], numobs)
            b = sigma.tilde.inv \%*\% mu.tilde + t(H.list[[l]][p1, 
                , ]) \%*\% (psi.list.inv[[l]][p1, , ]) \%*\% (t(yy) - 
                matrix(mu.list[[l]][, p1], r[l], numobs))
            chsi = ginv(A)
            if (!isSymmetric(chsi)) 
                chsi = makeSymm(chsi)
            roy = chsi \%*\% b
            roy.quadro = array(apply(roy, 2, function(x) x \%*\% 
                t(x)), c(r[l + 1], r[l + 1], numobs))
            zz2[, , , p1, p2, p3] = aperm(array(chsi, c(r[l + 
                1], r[l + 1], numobs)) + roy.quadro, c(3, 1, 
                2))
            z2.one[, , p1, p2, p3] = rmvnorm(numobs, rep(0, r[l + 
                1]), chsi) + t(roy)
            z2[, , p1, p2, p3] = t(roy)
        }
        for (i1 in 1:k[l + 1]) for (i2 in 1:k[l + 2]) {
            prob = ps.y.list[[l + 1]][, i1, drop = FALSE] * ps.y.list[[l + 
                2]][, i2, drop = FALSE]
            z = z + array(z2[, , , i1, i2, drop = FALSE] * array(rowSums(prob), 
                c(numobs, r[l + 1], k[l], 1, 1)), c(numobs, r[l + 
                1], k[l]))
            z.one = z.one + array(z2.one[, , , i1, i2, drop = FALSE] * 
                array(rowSums(prob), c(numobs, r[l + 1], k[l], 
                  1, 1)), c(numobs, r[l + 1], k[l]))
            zz = zz + array(zz2[, , , , i1, i2, drop = FALSE] * 
                array(rowSums(prob), c(numobs, r[l + 1], r[l + 
                  1], k[l], 1, 1)), c(numobs, r[l + 1], r[l + 
                1], k[l]))
        }
        z.list[[l]] = aperm(z.one, c(3, 1, 2))
        out = compute.est(k[l], r[l], r[l + 1], ps.y.list[[l]], 
            yy, aperm(z, c(3, 1, 2)), aperm(zz, c(4, 2, 3, 1)), 
            mu.list[[l]])
        H.list[[l]] = out$H
        psi.list[[l]] = out$psi
        psi.list.inv[[l]] = out$psi.inv
        mu.list[[l]] = out$mu
        w.list[[l]] = out$w
        l = 2
        yy = matrix(0, numobs, r[l])
        zz = z.list[[l - 1]]
        for (i in 1:k[l - 1]) {
            yy = yy + matrix(zz[i, , , drop = FALSE] * array(rowSums(ps.y.list[[l - 
                1]][, i, drop = FALSE]), c(1, numobs, r[l])), 
                numobs, r[l])
        }
        z2 <- z2.one <- array(0, c(numobs, r[l + 1], k[l], k[l + 
            1]))
        zz2 <- array(0, c(numobs, r[l + 1], r[l + 1], k[l], k[l + 
            1]))
        z <- z.one <- array(0, c(numobs, r[l + 1], k[l]))
        zz <- array(0, c(numobs, r[l + 1], r[l + 1], k[l]))
        for (p1 in 1:k[l]) for (p2 in 1:k[l + 1]) {
            A = ginv(H.list[[l + 1]][p2, , ] \%*\% t(H.list[[l + 
                1]][p2, , ]) + psi.list[[l + 1]][p2, , ]) + t(H.list[[l]][p1, 
                , ]) \%*\% (psi.list.inv[[l]][p1, , ]) \%*\% (H.list[[l]][p1, 
                , ])
            b = ginv(H.list[[l + 1]][p2, , ] \%*\% t(H.list[[l + 
                1]][p2, , ]) + psi.list[[l + 1]][p2, , ]) \%*\% 
                matrix(mu.list[[l + 1]][, p2], r[l + 1], numobs) + 
                t(H.list[[l]][p1, , ]) \%*\% (psi.list.inv[[l]][p1, 
                  , ]) \%*\% (t(yy) - matrix(mu.list[[l]][, p1], 
                  r[l], numobs))
            chsi = ginv(A)
            if (!isSymmetric(chsi)) 
                chsi = makeSymm(chsi)
            roy = chsi \%*\% b
            roy.quadro = array(apply(roy, 2, function(x) x \%*\% 
                t(x)), c(r[l + 1], r[l + 1], numobs))
            zz2[, , , p1, p2] = aperm(array(chsi, c(r[l + 1], 
                r[l + 1], numobs)) + roy.quadro, c(3, 1, 2))
            z2.one[, , p1, p2] = rmvnorm(numobs, rep(0, r[l + 
                1]), chsi) + t(roy)
            z2[, , p1, p2] = t(roy)
        }
        for (i in 1:k[l + 1]) {
            prob = ps.y.list[[l + 1]][, i, drop = FALSE]
            z = z + array(z2[, , , i, drop = FALSE] * array(rowSums(prob), 
                c(numobs, r[l + 1], k[l], 1)), c(numobs, r[l + 
                1], k[l]))
            z.one = z.one + array(z2.one[, , , i, drop = FALSE] * 
                array(rowSums(prob), c(numobs, r[l + 1], k[l], 
                  1)), c(numobs, r[l + 1], k[l]))
            zz = zz + array(zz2[, , , , i, drop = FALSE] * array(rowSums(prob), 
                c(numobs, r[l + 1], r[l + 1], k[l], 1)), c(numobs, 
                r[l + 1], r[l + 1], k[l]))
        }
        z.list[[l]] = aperm(z.one, c(3, 1, 2))
        out = compute.est(k[l], r[l], r[l + 1], ps.y.list[[l]], 
            yy, aperm(z, c(3, 1, 2)), aperm(zz, c(4, 2, 3, 1)), 
            mu.list[[l]])
        H.list[[l]] = out$H
        psi.list[[l]] = out$psi
        psi.list.inv[[l]] = out$psi.inv
        mu.list[[l]] = out$mu
        w.list[[l]] = out$w
        l = 3
        yy = matrix(0, numobs, r[l])
        zz = z.list[[l - 1]]
        for (i in 1:k[l - 1]) {
            yy = yy + matrix(zz[i, , , drop = FALSE] * array(rowSums(ps.y.list[[l - 
                1]][, i, drop = FALSE]), c(1, numobs, r[l])), 
                numobs, r[l])
        }
        z <- z.one <- array(0, c(numobs, r[l + 1], k[l]))
        zz <- array(0, c(numobs, r[l + 1], r[l + 1], k[l]))
        for (p1 in 1:k[l]) {
            A = diag(r[l + 1]) + t(H.list[[l]][p1, , ]) \%*\% (psi.list.inv[[l]][p1, 
                , ]) \%*\% (H.list[[l]][p1, , ])
            b = t(H.list[[l]][p1, , ]) \%*\% (psi.list.inv[[l]][p1, 
                , ]) \%*\% (t(yy) - matrix(mu.list[[l]][, p1], 
                r[l], numobs))
            chsi = ginv(A)
            if (!isSymmetric(chsi)) 
                chsi = makeSymm(chsi)
            roy = chsi \%*\% b
            roy.quadro = array(apply(roy, 2, function(x) x \%*\% 
                t(x)), c(r[l + 1], r[l + 1], numobs))
            zz[, , , p1] = aperm(array(chsi, c(r[l + 1], r[l + 
                1], numobs)) + roy.quadro, c(3, 1, 2))
            z.one[, , p1] = t(roy) + rmvnorm(numobs, rep(0, r[l + 
                1]), chsi)
            z[, , p1] = t(roy)
        }
        z.list[[l]] = aperm(z.one, c(3, 1, 2))
        out = compute.est(k[l], r[l], r[l + 1], ps.y.list[[l]], 
            yy, aperm(z, c(3, 1, 2)), aperm(zz, c(4, 2, 3, 1)), 
            mu.list[[l]])
        H.list[[l]] = out$H
        psi.list[[l]] = out$psi
        psi.list.inv[[l]] = out$psi.inv
        mu.list[[l]] = out$mu
        w.list[[l]] = out$w
        out = compute.lik(y, numobs, k, mu.list, H.list, psi.list, 
            w.list)
        py = out$py
        ps.y = out$ps.y
        ps.y.list = out$ps.y.list
        k.comb = out$k.comb
        s = out$s
        lik <- sum(log(py))
        likelihood <- c(likelihood, lik)
        if (hh < 5) 
            ratio <- 2 * eps
        if (hh > 5) 
            ratio <- (ma(likelihood, 5)[hh + 1] - ma(likelihood, 
                5)[hh])/abs(ma(likelihood, 5)[hh])
    }
    h = 0
    for (j in 1:layers) h = h + (k[j] - 1) + (r[j] * r[j + 1]) * 
        k[j] + r[j] * k[j] + k[j] * r[j]
    if (layers > 1) 
        for (j in 2:layers) h = h - (r[j] * k[j] * (r[j] - 1)/2)
    bic <- -2 * lik + h * log(numobs)
    aic <- -2 * lik + 2 * h
    EN <- entr(ps.y.list[[1]])
    clc <- -2 * lik + 2 * EN
    icl.bic <- -2 * lik + 2 * EN + h * log(numobs)
    out <- list(H = H.list, w = w.list, mu = mu.list, psi = psi.list, 
        likelihood = likelihood, bic = bic, aic = aic, clc = clc, 
        s = s, icl.bic = icl.bic, h = h, ps.y, ps.y)
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
