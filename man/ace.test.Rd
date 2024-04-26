\name{ace.test}
\alias{ace.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ACE permutation test of independence
%%  ~~function to do ... ~~
}
\description{
Performs a permutation test of independence using ace in package acepack.
ace stands for alternating conditional expectations.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ace.test(x, y = NULL, nperm = 999)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a numeric vector, or a matrix or data frame with two columns.
%%     ~~Describe \code{x} here~~
}
  \item{y}{NULL (default) or a vector with same length as x.
%%     ~~Describe \code{y} here~~
}
  \item{nperm}{number of permutations.
%%     ~~Describe \code{nperm} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{ace }{the value of the test statistic.}
  \item{pval }{the p-value of the test.}
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
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{lcor.test}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
n <- 200
x <- matrix(rnorm(n*2), n)
nu <- 2
y <- x / sqrt(rchisq(n, nu)/nu) #multivariate t
cor.test(y[,1], y[,2], method = "spearman")
ace.test(y)
}
