\name{compute.lik}
\alias{compute.lik}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
compute.lik(y, numobs, k, mu.list, H.list, psi.list, w.list)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{numobs}{
%%     ~~Describe \code{numobs} here~~
}
  \item{k}{
%%     ~~Describe \code{k} here~~
}
  \item{mu.list}{
%%     ~~Describe \code{mu.list} here~~
}
  \item{H.list}{
%%     ~~Describe \code{H.list} here~~
}
  \item{psi.list}{
%%     ~~Describe \code{psi.list} here~~
}
  \item{w.list}{
%%     ~~Describe \code{w.list} here~~
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
function (y, numobs, k, mu.list, H.list, psi.list, w.list) 
{
    value = 1e-20
    layers = length(k)
    p = ncol(y)
    py <- matrix(0, numobs)
    tot.k = prod(k)
    py.s <- matrix(0, numobs, tot.k)
    pys <- matrix(0, numobs, tot.k)
    k.comb = apply(t(k), 2, function(x) 1:x)
    if (is.list(k.comb)) 
        k.comb = expand.grid(k.comb)
    if (is.matrix(k.comb)) 
        k.comb = expand.grid(split(t(k.comb), 1:ncol(k.comb)))
    if (prod(k) == 1) 
        k.comb = matrix(k.comb, nrow = 1)
    for (i in 1:tot.k) {
        mu.tot = mu.list[[1]][, k.comb[i, 1]]
        var.tot = psi.list[[1]][k.comb[i, 1], , ]
        w.tot = w.list[[1]][k.comb[i, 1]]
        if (layers > 1) 
            for (l in 2:layers) {
                tot.H = diag(p)
                for (m in 1:(l - 1)) tot.H = tot.H \%*\% H.list[[m]][k.comb[i, 
                  m], , ]
                mu.tot = mu.tot + tot.H \%*\% mu.list[[l]][, k.comb[i, 
                  l]]
                var.tot = var.tot + tot.H \%*\% (H.list[[l]][k.comb[i, 
                  l], , ] \%*\% t(H.list[[l]][k.comb[i, l], , ]) + 
                  psi.list[[l]][k.comb[i, l], , ]) \%*\% t(tot.H)
                w.tot = w.tot * w.list[[l]][k.comb[i, l]]
            }
        if (!is.positive.definite(var.tot)) 
            var.tot = make.positive.definite(var.tot)
        py.s[, i] <- dmvnorm(y, c(mu.tot), as.matrix(var.tot), 
            log = TRUE)
        if (w.tot == 0) 
            w.tot = 10^(-320)
        pys[, i] = log(w.tot) + py.s[, i]
    }
    cc = 705 - max(pys)
    pys = pys + cc
    py = rowSums(exp(pys))
    ps.y = exp(pys)/matrix(py, numobs, tot.k)
    ps.y = ifelse(is.na(ps.y), 1/k, ps.y)
    py = exp(-cc) * py
    s = matrix(0, numobs, layers)
    ps.y.list = NULL
    for (l in 1:layers) {
        psi.y = matrix(0, numobs, k[l])
        for (i in 1:k[l]) {
            index = (k.comb[, l] == i)
            psi.y[, i] = rowSums(ps.y[, index, drop = FALSE])
        }
        s[, l] = apply(psi.y, 1, which.max)
        ps.y.list[[l]] = psi.y
    }
    return(list(py = py, py.s = py.s, ps.y = ps.y, k.comb = k.comb, 
        s = s, ps.y.list = ps.y.list))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
