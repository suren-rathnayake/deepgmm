\name{model.selection}
\alias{model.selection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
model.selection(y, k, layers, seeds = 3, it = 50, eps = 0.001, criterion = "BIC")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{k}{
%%     ~~Describe \code{k} here~~
}
  \item{layers}{
%%     ~~Describe \code{layers} here~~
}
  \item{seeds}{
%%     ~~Describe \code{seeds} here~~
}
  \item{it}{
%%     ~~Describe \code{it} here~~
}
  \item{eps}{
%%     ~~Describe \code{eps} here~~
}
  \item{criterion}{
%%     ~~Describe \code{criterion} here~~
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
function (y, k, layers, seeds = 3, it = 50, eps = 0.001, criterion = "BIC") 
{
    bic.best = Inf
    aic.best = Inf
    p = dim(y)[2]
    pp = round(p/2, 0)
    ppp = round(p/3, 0)
    pppp = round(p/4, 0)
    if (layers == 1) {
        r = c(1:pp)
        bic = array(NA, c(seeds, pp))
        bic.best = Inf
        aic = array(NA, c(seeds, pp))
        aic.best = Inf
        for (i in 1:seeds) for (rr in 1:pp) {
            out = try(deep.sem(y, 1, k, rr, it = 50, eps = 0.001, 
                seed = i))
            if (!is.character(out)) {
                if (criterion == "BIC") 
                  if (out$bic < bic.best) {
                    out.best = out
                    bic.best = out$bic
                  }
                if (criterion == "AIC") 
                  if (out$aic < aic.best) {
                    out.best = out
                    aic.best = out$aic
                  }
                bic[i, rr] = out$bic
                aic[i, rr] = out$aic
            }
        }
        if (criterion == "BIC") 
            index = which(bic == min(bic), arr.ind = TRUE)[1, 
                ]
        if (criterion == "AIC") 
            index = which(aic == min(aic), arr.ind = TRUE)[1, 
                ]
        message("Best Fit: \n")
        cat(paste("Seed=", index[1], " r=", index[2], " BIC:", 
            round(out.best$bic, 2), " AIC:", round(out.best$aic, 
                2)))
    }
    if (layers == 2) {
        r = as.matrix(expand.grid(1:pp, 1:ppp))
        r = r[(r[, 1]) > (r[, 2]), ]
        k = rbind(c(k, 1), c(k, 2), c(k, 3))
        bic = array(NA, c(seeds, nrow(k), nrow(r)))
        bic.best = Inf
        aic = array(NA, c(seeds, nrow(k), nrow(r)))
        aic.best = Inf
        for (i in 1:seeds) for (kk in 1:nrow(k)) for (rr in 1:nrow(r)) {
            out = try(deep.sem(y, 2, k[kk, ], r[rr, ], it = 50, 
                eps = 0.001, seed = i))
            if (!is.character(out)) {
                if (criterion == "BIC") 
                  if (out$bic < bic.best) {
                    out.best = out
                    bic.best = out$bic
                  }
                if (criterion == "AIC") 
                  if (out$aic < aic.best) {
                    out.best = out
                    aic.best = out$aic
                  }
                bic[i, kk, rr] = out$bic
                aic[i, kk, rr] = out$aic
            }
        }
        if (criterion == "BIC") 
            index = which(bic == min(bic), arr.ind = TRUE)[1, 
                ]
        if (criterion == "AIC") 
            index = which(aic == min(aic), arr.ind = TRUE)[1, 
                ]
        message("Best Fit: \n")
        cat(paste("Seed=", index[1], " k=", paste(k[index[2], 
            ], collapse = " "), " r=", paste(r[index[3], ], collapse = " "), 
            " BIC:", round(out.best$bic, 2), " AIC:", round(out.best$aic, 
                2)))
    }
    if (layers == 3) {
        r = as.matrix(expand.grid(1:pp, 1:ppp, 1:ppp))
        r = r[((r[, 1]) > (r[, 2])) & ((r[, 2]) > (r[, 3])), 
            ]
        k = rbind(c(k, 1, 1), c(k, 2, 1), c(k, 3, 1), c(k, 1, 
            2), c(k, 2, 2), c(k, 3, 2), c(k, 1, 3), c(k, 2, 3), 
            c(k, 3, 3))
        bic = array(NA, c(seeds, nrow(k), nrow(r)))
        bic.best = Inf
        aic = array(NA, c(seeds, nrow(k), nrow(r)))
        aic.best = Inf
        for (i in 1:seeds) for (kk in 1:nrow(k)) for (rr in 1:nrow(r)) {
            out = try(deep.sem(y, 3, k[kk, ], r[rr, ], it = 50, 
                eps = 0.001, seed = i))
            if (!is.character(out)) {
                if (criterion == "BIC") 
                  if (out$bic < bic.best) {
                    out.best = out
                    bic.best = out$bic
                  }
                if (criterion == "AIC") 
                  if (out$aic < aic.best) {
                    out.best = out
                    aic.best = out$aic
                  }
                bic[i, kk, rr] = out$bic
                aic[i, kk, rr] = out$aic
            }
        }
        if (criterion == "BIC") 
            index = which(bic == min(bic), arr.ind = TRUE)[1, 
                ]
        if (criterion == "AIC") 
            index = which(aic == min(aic), arr.ind = TRUE)[1, 
                ]
        message("Best Fit: \n")
        cat(paste("Seed=", index[1], " k=", paste(k[index[2], 
            ], collapse = " "), " r=", paste(r[index[3], ], collapse = " "), 
            " BIC:", round(out.best$bic, 2), " AIC:", round(out.best$aic, 
                2)))
    }
    out = list(fit = out.best, bic = bic, aic = aic)
    invisible(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
