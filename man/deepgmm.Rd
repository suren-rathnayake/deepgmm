\name{deepgmm}
\alias{deepgmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Fits Deep Gaussian Mixture Models Using
 Stochastic EM algorithm.
}
\description{
  Fits deep Gaussian mixture model to a multivariate data.
}
\usage{
deepgmm(y, layers, k, r = rep(1, layers),
        it = 50, eps = 0.001, init = "kmeans", method = "factanal")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
   A matrix or a data frame of which rows correspond to
   observations and columns to variables.
}
  \item{layers}{
   The number of layers in the deep Gaussian mixture model.
   Admitted values are 1, 2 or 3.
}
  \item{k}{
   A vector of integers of length \code{layers}
   containing the number of groups in the different layers.
}
  \item{r}{
   A vector of integers of length \code{layers}
   containing the dimensions at the different layers.
   Dimension of the layers must be in decreasing
   size.
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
\item{method}{
 To determine how the initial parameter values are computed. See Details.
}
}
\details{
A Deep Gaussian Mixture model is a network of multiple
layers of latent variables, where, at each layer, the variables follow a
mixture of Gaussian distributions.
Thus, the deep mixture model consists of a set of nested
mixtures of linear models, which globally provide a nonlinear model
able to describe the data in a very  flexible way. In order to avoid over
parameterized solutions, dimension reduction by factor models can be
applied at each layer of the architecture thus resulting in deep mixtures
of factor analysers.
}
\value{
An object of class \code{"dgmm"} containing fitted values.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
    Viroli, C. and McLachlan, G.J. (2018).
    Deep Gaussian mixture models. Statistics and Computing.
    (Advance Access, published 1 December, 2017). To appear.
    Preprint arXiv:1711.06929.
}
\author{
 Cinzia Viroli, Geoffrey J. McLachlan
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}
\examples{
\dontrun{
library(gclus)
data(wine)

# Scale variables
y <- scale(wine[, -1])
cls <- wine[, 1]

## fit a DGMM with two layers
layers <- 2

## number of groups in the different layers
k <- c(3, 2)
# 3 is the number of clusters at the observed level,
# 2 is the number of clusters at the latent level

## dimensions at the different layers
r <- c(5, 1)

set.seed(1)
fit <- deepgmm(y, layers, k, r, it=250, eps=0.001)
fit

summary(fit)
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}% use one of  RShowDoc("KEYWORDS")
\keyword{models}% __ONLY ONE__ keyword per line
\keyword{multivariate}
