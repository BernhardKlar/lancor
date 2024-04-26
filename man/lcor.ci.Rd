\name{lcor.ci}
\alias{lcor.ci}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{confidence intervals for the Lancaster correlation coefficient
%%  ~~function to do ... ~~
}
\description{
Computes confidence intervals for the Lancaster correlation coefficient.
Lancaster correlation is a bivariate measures of dependence.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
lcor.ci(x, y = NULL, conf.level = 0.95, type = c("rank", "linear"), con = TRUE,
        R = 1000, method = c("plugin", "boot", "pretest"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a numeric vector, or a matrix or data frame with two columns.
%%     ~~Describe \code{x} here~~
}
  \item{y}{NULL (default) or a vector with same length as x.
%%     ~~Describe \code{y} here~~
}
  \item{conf.level}{
%%     ~~Describe \code{conf.level} here~~
}
  \item{type}{a character string indicating which lancaster correlation is to be computed. One of "rank" (default), or "linear": can be abbreviated.
%%     ~~Describe \code{type} here~~
}
  \item{con}{logical; if TRUE (default), conservative asymptotic confidence intervals are computed.
%%     ~~Describe \code{con} here~~
}
  \item{R}{number of bootstrap replications.
%%     ~~Describe \code{R} here~~
}
  \item{method}{a character string indicating how the asymptotic covariance matrix is computed if type ="linear". One of "plugin" (default), "boot" or "symmetric": can be abbreviated.
%%     ~~Describe \code{method} here~~
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
\verb{lcor.ci} returns a vector containing the left and right endpoint of the confidence interval.
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

\seealso{\code{\link{lcor}}, \code{\link{lcor.comp}}, \code{\link{lcor.test}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
n <- 1000
x <- matrix(rnorm(n*2), n)
nu <- 2
y <- x / sqrt(rchisq(n, nu)/nu) # multivariate t
lcor(y, type = "rank")
lcor.ci(y, type = "rank")
}
