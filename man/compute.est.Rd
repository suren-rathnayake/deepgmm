\name{compute.est}
\alias{compute.est}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
compute.est(k, p, qq, ps.y, y, Ez, Ezz, mu)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{k}{
%%     ~~Describe \code{k} here~~
}
  \item{p}{
%%     ~~Describe \code{p} here~~
}
  \item{qq}{
%%     ~~Describe \code{qq} here~~
}
  \item{ps.y}{
%%     ~~Describe \code{ps.y} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{Ez}{
%%     ~~Describe \code{Ez} here~~
}
  \item{Ezz}{
%%     ~~Describe \code{Ezz} here~~
}
  \item{mu}{
%%     ~~Describe \code{mu} here~~
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
function (k, p, qq, ps.y, y, Ez, Ezz, mu) 
{
    value = 1e-20
    H = array(0, c(k, p, qq))
    psi = psi.inv = array(0, c(k, p, p))
    numobs = dim(y)[1]
    w = 0
    for (i in 1:k) {
        if (qq > 1) 
            EEzz = apply((Ezz[i, , , ] * (aperm(array(t(t(ps.y[, 
                i])), c(numobs, qq, qq)), c(2, 3, 1)))), 1, rowSums)/sum(ps.y[, 
                i])
        if (qq == 1) 
            EEzz = Ezz[i, , , ] \%*\% (ps.y[, i])/sum(ps.y[, i])
        H[i, , ] = (t((y - t(matrix(mu[, i], p, numobs))) * matrix(ps.y[, 
            i], numobs, p)) \%*\% (matrix(Ez[i, , ], nrow = numobs) * 
            matrix(ps.y[, i], numobs, qq))) \%*\% ginv(EEzz)/sum(ps.y[, 
            i])
        H[i, , ] = ifelse(is.na(H[i, , ]), mean(H[i, , ], na.rm = TRUE), 
            H[i, , ])
        psi[i, , ] = (t((y - t(matrix(mu[, i], p, numobs))) * 
            matrix(ps.y[, i], numobs, p)) \%*\% (matrix((y - t(matrix(mu[, 
            i], p, numobs))) * matrix(ps.y[, i], numobs, p), 
            ncol = p)) - t((y - t(matrix(mu[, i], p, numobs))) * 
            matrix(ps.y[, i], numobs, p)) \%*\% (matrix(Ez[i, , 
            ], nrow = numobs) * matrix(ps.y[, i], numobs, qq)) \%*\% 
            t(H[i, , ]))
        psi[i, , ] = psi[i, , , drop = FALSE]/sum(ps.y[, i])
        psi[i, , ] = ifelse((psi[i, , ] == 0), value, psi[i, 
            , ])
        if (p > 1) {
            psi[i, , ] <- diag(diag(psi[i, , ]))
            psi.inv[i, , ] <- diag(1/diag(psi[i, , ]))
        }
        psi[i, , ] = ifelse(is.na(psi[i, , ]), 1, psi[i, , ])
        psi.inv[i, , ] = ifelse(is.na(psi.inv[i, , ]), 1, psi.inv[i, 
            , ])
        mu[, i] = colSums(matrix(ps.y[, i], numobs, p) * (y - 
            t(matrix(H[i, , ], ncol = qq) \%*\% t(Ez[i, , ]))))/sum(ps.y[, 
            i])
        w[i] = mean(ps.y[, i])
    }
    return(list(H = H, mu = mu, psi = psi, psi.inv = psi.inv, 
        w = w))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
