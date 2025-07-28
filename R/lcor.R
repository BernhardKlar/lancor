#' Lancaster correlation
#' 
#' @description 
#' Computes the Lancaster correlation coefficient.
#' 
#' @param x a numeric vector, or a matrix or data frame with two columns.
#' @param y NULL (default) or a vector with same length as x.
#' @param type a character string indicating which lancaster correlation is to be computed. One of "rank" (default), or "linear": can be abbreviated.
#'
#' @return 
#' the sample Lancaster correlation.
#' 
#' @author Hajo Holzmann, Bernhard Klar
#' 
#' @references
#' Holzmann, Klar (2024). "Lancester correlation - a new dependence measure linked to maximum correlation". \url{https://arxiv.org/abs/2303.17872}
#' 
#' @seealso \code{\link{lcor.comp}, \link{lcor.ci}, \link{lcor.test}}
#'
#' @examples 
#' Sigma <- matrix(c(1,0.1,0.1,1), ncol=2)
#' R <- chol(Sigma)
#' n <- 1000
#' x <- matrix(rnorm(n*2), n)
#' lcor(x, type = "rank")
#' lcor(x, type = "linear")
#' 
#' x <- matrix(rnorm(n*2), n)
#' nu <- 2
#' y <- x / sqrt(rchisq(n, nu)/nu)
#' cor(y[,1], y[,2], method = "spearman")
#' lcor(y, type = "rank")
#'
#' @export 
lcor = function(x, y = NULL, type = c("rank", "linear")) {
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
    y = qnorm( (rank(y)-0.5) / length(x) )
    lc = max( abs(cor(x,y)), abs(cor(x^2,y^2)) )
    return(lc)
  }
  if (type == "linear") { 
    x = scale(x)
    y = scale(y)
    lc = max( abs(cor(x,y)), abs( cor(x^2,y^2) ) )
    return(lc)
  }
}
