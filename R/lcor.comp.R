#' Lancaster correlation and its components
#' 
#' @description 
#' Computes the Lancaster correlation coefficient and its components.
#'
#' @param x a numeric vector, or a matrix or data frame with two columns.
#' @param y NULL (default) or a vector with same length as x. 
#' @param type a character string indicating which lancaster correlation is to be computed. One of "rank" (default), or "linear": can be abbreviated.
#' @param plot logical; if TRUE, scatterplots of the transformed x and y values and of their squares are drawn.
#'
#' @return 
#' a vector containing the two components rho1 and rho2 and the sample Lancaster correlation.
#' 
#' @author Hajo Holzmann, Bernhard Klar
#' 
#' @references
#' Holzmann, Klar (2024). "Lancester correlation - a new dependence measure linked to maximum correlation". \url{https://arxiv.org/abs/2303.17872}
#' 
#' @seealso \code{\link{lcor}, \link{lcor.comp}, \link{lcor.test}}
#'
#' @examples
#' Sigma <- matrix(c(1,0.1,0.1,1), ncol=2)
#' R <- chol(Sigma)
#' n <- 1000
#' x <- matrix(rnorm(n*2), n)
#' nu <- 8
#' y <- x / sqrt(rchisq(n, nu)/nu) #multivariate t
#' cor(y[,1], y[,2])
#' lcor.comp(y, type = "linear")
#' 
#' x <- matrix(rnorm(n*2), n)
#' nu <- 2
#' y <- x / sqrt(rchisq(n, nu)/nu) #multivariate t
#' cor(y[,1], y[,2], method = "spearman")
#' lcor.comp(y, type = "rank", plot = TRUE)
#' 
#'
#' @export 
lcor.comp = function(x, y = NULL, type = c("rank", "linear"), plot = FALSE) {
  if (is.data.frame(x)) x = as.matrix(x)
  if (!is.matrix(x) && is.null(y))
    stop("supply both 'x' and 'y' or a matrix-like 'x'")
  if (is.matrix(x)) {
    y = x[,2]
    x = x[,1]
    }
  type = match.arg(type)
  if (!(type %in% c("rank", "linear")))
    stop("type for lcor can only be rank or linear")
  if (type == "rank") {
    x = qnorm( (rank(x)-0.5) / length(x) )
    y = qnorm( (rank(y)-0.5) / length(y) )
    rho1 = cor(x,y)
    rho2 = cor(x^2,y^2)
    lc = max( abs(rho1), abs(rho2) )
  }
  if (type == "linear") {
    x = scale(x)
    y = scale(y)
    rho1 = cor(x,y)
    rho2 = cor(x^2,y^2)
    lc = max( abs(rho1), abs(rho2) )
  }
  #
  if (plot==TRUE) {
    oldpar = par(no.readonly = TRUE)
    on.exit(par(oldpar))
    par(mfrow=c(1,2))
    plot(y ~ x, pch=16, xlab="x.standardized", ylab="y.standardized")
    plot(y^2 ~ I(x^2), cex=0.7, pch=16, xlab="square of x.standardized",
         ylab="square of y.standardized")
  }
  return( setNames( c(rho1, rho2, lc), c("rho1", "rho2", "lcor") ) )
}
