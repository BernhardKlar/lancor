\name{lcor}
\alias{lcor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Lancaster correlation
%%  ~~function to do ... ~~
}
\description{
Computes the Lancaster correlation coefficient.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
lcor(x, y = NULL, type = c("rank", "linear"))
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
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\verb{lcor} returns the sample Lancaster correlation.
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

\seealso{\code{\link{lcor.comp, lcor.ci, lcor.test}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
Sigma <- matrix(c(1,0.1,0.1,1), ncol=2)
R <- chol(Sigma)
n <- 1000
x <- matrix(rnorm(n*2), n) %*% R #multivariate normal
lcor(x, type = "rank")
lcor(x, type = "linear")

x <- matrix(rnorm(n*2), n)
nu <- 2
y <- x / sqrt(rchisq(n, nu)/nu)
cor(y[,1], y[,2], method = "spearman")
lcor(y, type = "rank")
}
