\name{model_selection}
\alias{model_selection}

\title{
 Function to compare different models
}
\description{
  Compares different models and return the best one selected according to
  criterion (BIC or AIC).
}
\usage{
model_selection(y, layers, g, seeds = 3, it = 50, eps = 0.001,
        init = "kmeans", init_est = "factanal", criterion = "BIC")
}
\arguments{
  \item{y}{
   A matrix or a data frame in which rows correspond to
   observations and columns to variables.
}
  \item{layers}{
   The number of layers in the deep Gaussian mixture model.
   Admitted values are 1, 2 or 3.
}
  \item{g}{
  The number of clusters.
}
  \item{seeds}{
   Numeric vector containing seeds to try.
}
  \item{it}{
   Maximum number of EM iterations.
}
  \item{eps}{
   The EM algorithm terminates the relative increment of the log-likelihod
   falls below this value.
}
  \item{init}{
   Initial paritioning of the observations to determine initial
   parameter values. See Details.
}
\item{init_est}{
 To determine how the initial parameter values are computed. See Details.
}
\item{criterion}{
 Model selection criterion, either \code{"AIC"} of \code{"BIC"}.
}
}
\details{
Compares different models and return the best one selected according to
criterion (BIC or AIC). One can use diffefrent number of seeds.
}
\value{
A list containing
an object of class \code{"dgmm"} containing fitted values
and list of BIC and AIC values.
}
\references{
    Viroli, C. and McLachlan, G.J. (2019). Deep Gaussian mixture models. 
    Statistics and Computing 29, 43-51.
}
% \author{
% %%  ~~who you are~~
% }
% \note{
% %%  ~~further notes~~
% }
\examples{
\donttest{
layers <- 2
k <- c(3, 4) 
r <- c(3, 2)
it <- 50
eps <- 0.001
y <- scale(mtcars)

sel <- model_selection(y, layers, 3, seeds = 1, it = 250, eps = 0.001)
sel

summary(sel)
}
}

\keyword{cluster}% use one of  RShowDoc("KEYWORDS")
\keyword{models}% __ONLY ONE__ keyword per line
\keyword{multivariate}
