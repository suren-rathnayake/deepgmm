\name{mapClass}
\alias{mapClass}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mapClass(a, b)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{
%%     ~~Describe \code{a} here~~
}
  \item{b}{
%%     ~~Describe \code{b} here~~
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
function (a, b) 
{
    l <- length(a)
    x <- y <- rep(NA, l)
    if (l != length(b)) {
        warning("unequal lengths")
        return(x)
    }
    aChar <- as.character(a)
    bChar <- as.character(b)
    Tab <- table(a, b)
    Ua <- dimnames(Tab)[[1]]
    Ub <- dimnames(Tab)[[2]]
    aTOb <- rep(list(Ub), length(Ua))
    names(aTOb) <- Ua
    bTOa <- rep(list(Ua), length(Ub))
    names(bTOa) <- Ub
    k <- nrow(Tab)
    Map <- rep(0, k)
    Max <- apply(Tab, 1, max)
    for (i in 1:k) {
        I <- match(Max[i], Tab[i, ], nomatch = 0)
        aTOb[[i]] <- Ub[I]
    }
    if (is.numeric(b)) 
        aTOb <- lapply(aTOb, as.numeric)
    k <- ncol(Tab)
    Map <- rep(0, k)
    Max <- apply(Tab, 2, max)
    for (j in (1:k)) {
        J <- match(Max[j], Tab[, j])
        bTOa[[j]] <- Ua[J]
    }
    if (is.numeric(a)) 
        bTOa <- lapply(bTOa, as.numeric)
    list(aTOb = aTOb, bTOa = bTOa)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
