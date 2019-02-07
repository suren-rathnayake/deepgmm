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
        it = 250, eps = 0.001, init = "kmeans", init_est = "factanal")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
   A matrix or a data frame in which the rows correspond to
   observations and the columns to variables.
}
  \item{layers}{
   The number of layers in the deep Gaussian mixture model.
   Limited to 1, 2 or 3.
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
\item{init_est}{
 Procedure for computing the initial parameter values for the given initial partition of the data. See Details.
}
}
\details{
The deep Gaussian mixture model is an hierarchical model organized
in a multilayered architecture where, at each layer, 
the variables follow a mixture of Gaussian distributions.
This set of nested mixtures of linear models provides a globally 
nonlinear model that can model the data in a very flexible way. 
In order to avoid overparameterized solutions, 
dimension reduction by factor models can be applied at each layer of 
the architecture, thus resulting in deep mixtures of factor analyzers.

The data \code{y} must be a matrix or a data frame containing
numerical values, with no missing values. The rows must correspond to
observations and the columns to variables.

Presently, the maximum number of layers \code{layers} implemented
is 3. 

The ith element of \code{k} contain number of groups in the ith layer. Thus
the length \code{k} must equal to \code{layers}. 

The parameter vector \code{r} contains the latent variable dimension of 
each layer. 
Variables at different layers have progressively decreasing dimension,
\eqn{r_1}, \eqn{r_2}, \dots, \eqn{r_h}, where \eqn{p > r_1 > r_2 > 
\dots > r_h \geq 1}.

The EM algorithm used by \code{dgmm} requires initialization. 
The initialization is done by partitioning the dataset, 
and then estimating the initial values for model parameters 
based on the partition. 
There are three options available in \code{dgmm} for the
initial partitioning of the data;
random partitioning, clustering using the \emph{k}-means algorithm
and using the agglomerative hierarchical clustering.
With the \code{init = "random"}
the partitioning is done randomly,
with \code{init = "kmeans"} the data are partitioned using the
\emph{k}-means algorithm of "Hartigan-Wong", 
while with \code{init = "hclass"} agglomerative 
hierarchical clustering is done to obtain initial partitioning. 
For this option the dissimilarity structure is obtained using the 
Euclidean distance measure.

After the initial partitioning has been chosen, initial values of
the parameters in the component analyzers need to be 
calculated. There are two options available in \code{init_est}.
The default option, \code{init_est = "factanal"} provides estimates of the
initial parameters based on factor analysis. If 
\code{init_est = "ppca"} then mixtures of 
probabilistic principal component analyzers are fitted within each layer
to estimate values for model parameters. 
}
\value{
An object of class \code{"dgmm"} containing fitted values.
It contains
 \item{H}{A list in which the \emph{i}th element is the loading 
 matrix for the \emph{i}th layer}
 \item{w}{A list containing mixing proportions for each layer. (i.e. 
 the element w[[i]][j] contain the mixing proportion of the \emph{j}th
 component in the \emph{i} layer.)}
 \item{mu}{A list of matrices containing components means in the columns.
  (i.e. the element mu[[i]][, j] contain the component mean of the \emph{j}th
 component in the \emph{i} layer.)}
 \item{psi}{A list of arrays which contain covariance matrices for 
 the random error components of each component (i.e. the element 
 psi[[i]][j, ,, ] contain the error covariance matrix 
 for the \emph{j}th component in the \emph{i} layer.)} 
\item{lik}{The log-likelihood after each EM iteration}
\item{bic}{The Bayesian information criterion for the model fit}
\item{acl}{The Akaike information criterion for the model fit}
\item{clc}{The Classification likelihood information criterion for the model fit}
\item{icl.bic}{The integrated classification criterion for the model fit}
\item{s}{Clustering of the observations}
}
\references{
    Viroli, C. and McLachlan, G.J. (2019). Deep Gaussian mixture models. Statistics and Computing 29, 43-51.
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
