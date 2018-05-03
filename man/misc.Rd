\name{misc}
\alias{misc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
misc(classification, truth)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{classification}{
%%     ~~Describe \code{classification} here~~
}
  \item{truth}{
%%     ~~Describe \code{truth} here~~
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
function (classification, truth) 
{
    q <- function(map, len, x) {
        x <- as.character(x)
        map <- lapply(map, as.character)
        y <- sapply(map, function(x) x[1])
        best <- y != x
        if (all(len) == 1) 
            return(best)
        errmin <- sum(as.numeric(best))
        z <- sapply(map, function(x) x[length(x)])
        mask <- len != 1
        counter <- rep(0, length(len))
        k <- sum(as.numeric(mask))
        j <- 0
        while (y != z) {
            i <- k - j
            m <- mask[i]
            counter[m] <- (counter[m]\%\%len[m]) + 1
            y[x == name(map)[m]] <- map[[m]][counter[m]]
            temp <- y != x
            err <- sum(as.numeric(temp))
            if (err < errmin) {
                errmin <- err
                best <- temp
            }
            j <- (j + 1)\%\%k
        }
        best
    }
    if (any(isNA <- is.na(classification))) {
        classification <- as.character(classification)
        nachar <- paste(unique(classification[!isNA]), collapse = "")
        classification[isNA] <- nachar
    }
    MAP <- mapClass(classification, truth)
    len <- sapply(MAP[[1]], length)
    if (all(len) == 1) {
        CtoT <- unlist(MAP[[1]])
        I <- match(as.character(classification), names(CtoT), 
            nomatch = 0)
        one <- CtoT[I] != truth
    }
    else {
        one <- q(MAP[[1]], len, truth)
    }
    len <- sapply(MAP[[2]], length)
    if (all(len) == 1) {
        TtoC <- unlist(MAP[[2]])
        I <- match(as.character(truth), names(TtoC), nomatch = 0)
        two <- TtoC[I] != classification
    }
    else {
        two <- q(MAP[[2]], len, classification)
    }
    err <- if (sum(as.numeric(one)) > sum(as.numeric(two))) 
        as.vector(one)
    else as.vector(two)
    bad <- seq(along = classification)[err]
    errorRate = length(bad)/length(truth)
    return(errorRate)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
