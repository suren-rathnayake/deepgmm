\name{deepgmm}
\alias{deepgmm}
\title{
 Fits Deep Gaussian Mixture Models Using
 Stochastic EM algorithm.
}
\description{
  Fits a deep Gaussian mixture model to multivariate data.
}
\usage{
deepgmm(y, layers, k, r,
        it = 250, eps = 0.001, init = "kmeans", method = "factanal")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
   A matrix or a data frame of which the rows correspond to
   observations and the columns to variables.
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
   size. See details.
}
  \item{it}{
   Maximum number of EM iterations.
}
  \item{eps}{
   The EM algorithm terminates if the relative increment of the log-likelihood
   falls below this value.
}
  \item{init}{
   Initial partitioning of the observations to determine initial
   parameter values. See Details.
}
\item{method}{
 To determine how the initial parameter values are computed. See Details.
}
}
\details{
Deep Gaussian mixture model is hierarchical inference method organized
in a multi-layered architecture, where, at each layer, 
the variables follow a mixture of Gaussian distributions.
This set of nested mixtures of linear models provide a globally 
nonlinear model that can to describe the data in a very flexible way. 
In order to avoid over parameterized solutions, 
dimension reduction by factor models can be applied at each layer of 
the architecture thus resulting in deep mixtures of factor analyzers.

The data \code{y} must be a matrix or a data frame containing
numerical values, with no missing values. The rows must correspond to
observations and the columns to variables.

Presently, the maximum number of layers \code{layers} implemented
is 3. 

The ith element of \code{k} contain number of groups in the ith layer. Thus
the length \code{k} must equal to \code{layers}. 

The parameter \code{r} contain the latent variable dimension of each layer. 
Variable at different layers have progressively decreasing dimension,
\eqn{r_1}, \eqn{r_2}, \dots, \eqn{r_h}, where \eqn{p > r_1 > r_2 > 
\dots > r_h \geq 1}.

The EM algorithm used by \code{dgmm} require initialization. 
The initialization is done by partitioning the dataset.
For this, three options are available in \code{dgmm};
random partitioning, clustering using the k-means algorithm
and using the hierarchical clustering.
 With the \code{init = "random"},
the partitioning is done randomly,
With \code{init = "kmeans"}, the data are partitioned using the
k-means algorithm of "Hartigan-Wong", 
while \code{init = "hclass"}
hierarchical clustering is done to obtain initial partitioning. 
For this the dissimilarity structure obtained using the 
Euclidean distance measure.

By default, \code{method = "factanal"}, estimation of initial parameter
values for the SEM algorithm is done using an approach based on 
the maximum-likelihood factor analysis. Otherwise, mixtures of 
probabilistic principal component analysis is done within each layer
to estimate initial model parameters. 
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
