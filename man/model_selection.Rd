\name{model_selection}
\alias{model_selection}

\title{
 Function to compare different models
}
\description{
  Compare different models and return the best one selected according to
  criterion (BIC or AIC).
}
\usage{
model_selection(y, layers, g, seeds = 3, it = 50, eps = 0.001,
        init = "kmeans", method = "factanal", criterion = "BIC")
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
\item{method}{
 To determine how the initial parameter values are computed. See Details.
}
\item{criterion}{
 Model selection criterion, either \code{"AIC"} of \code{"BIC"}.
}
}
\details{
Compare different models and return the best one selected according to
criterion (BIC or AIC). One can use diffefrent number of seeds.
}
\value{
A list containing
an object of class \code{"dgmm"} containing fitted values
and list of BIC and AIC values.
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
% \author{
% %%  ~~who you are~~
% }
% \note{
% %%  ~~further notes~~
% }
%
% %% ~Make other sections like Warning with \section{Warning }{....} ~
%
% \seealso{
% %% ~~objects to See Also as \code{\link{help}}, ~~~
% }
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
k <- 3
# 3 is the number of clusters at the observed level,
# 2 is the number of clusters at the latent level

sel <- model_selection(y, layers, k, seeds = 3, it = 250, eps = 0.001)
sel

summary(sel)
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}% use one of  RShowDoc("KEYWORDS")
\keyword{models}% __ONLY ONE__ keyword per line
\keyword{multivariate}
