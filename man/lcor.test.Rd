\name{lcor.test}
\alias{lcor.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Lancaster correlation test
%%  ~~function to do ... ~~
}
\description{
Lancaster correlation test of bivariate independence.
Lancaster correlation is a bivariate measures of dependence.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
lcor.test(x, y = NULL, type = c("rank", "linear"), nperm = 999,
          method = c("permutation", "asymptotic", "symmetric"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a numeric vector, or a matrix or data frame with two columns.
%%     ~~Describe \code{x} here~~
}
  \item{y}{NULL (default) or a vector with same length as x.
%%     ~~Describe \code{y} here~~
}
  \item{type}{a character string indicating which lancaster correlation is to be computed. One of "rank" (default), or "linear": can be abbreviated.
%%     ~~Describe \code{type} here~~
}
  \item{nperm}{number of permutations.
%%     ~~Describe \code{nperm} here~~
}
  \item{method}{a character string indicating how the p-value is computed if type ="linear". One of "permutation" (default), "asymptotic" or "symmetric": can be abbreviated.
%%     ~~Describe \code{method} here~~
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
  \item{lcor }{the value of the test statistic.}
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
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{lcor}}, \code{\link{lcor.comp}}, \code{\link{lcor.ci}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
n <- 200
x <- matrix(rnorm(n*2), n)
nu <- 2
y <- x / sqrt(rchisq(n, nu)/nu)
cor.test(y[,1], y[,2], method = "spearman")
lcor.test(y, type = "rank")
}
