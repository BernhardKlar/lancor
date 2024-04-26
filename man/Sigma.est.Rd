\name{Sigma.est}
\alias{Sigma.est}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Covariance matrix of components of Lancaster correlation coefficient.
%%  ~~function to do ... ~~
}
\description{
Estimate of covariance matrix of the two components of Lancaster correlation.
Lancaster correlation is a bivariate measures of dependence.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Sigma.est(xx)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xx}{a matrix or data frame with two columns.
%%     ~~Describe \code{xx} here~~
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
\verb{Sigma.est} returns the estimated covariance matrix.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Holzmann, Klar (2024) Lancester correlation - a new dependence
measure linked to maximum correlation. arXiv:2303.17872
%% ~put references to the literature/web site here ~
}
\author{Hajo Holzmann, Bernhard Klar
%%  ~~who you are~~
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{lcor.ci}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
Sigma <- matrix(c(1,0.1,0.1,1), ncol=2)
R <- chol(Sigma)
n <- 1000
x <- matrix(rnorm(n*2), n) %*% R
nu <- 8
y <- x / sqrt(rchisq(n, nu)/nu) #multivariate t
Sigma.est(y)
}
