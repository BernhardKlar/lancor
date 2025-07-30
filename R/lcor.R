#' Lancaster correlation
#' 
#' @description 
#' Computes the Lancaster correlation coefficient.
#' 
#' @details
#' Let \eqn{F_X} and \eqn{F_Y} be the distribution functions of \eqn{X} and \eqn{Y}, and define
#' \deqn{X^* = \Phi^{-1}(F_X(X)), \quad Y^* = \Phi^{-1}(F_Y(Y)),}
#' where \eqn{\Phi^{-1}} is the standard normal quantile function. Then
#' \deqn{\rho_L(X,Y) = \max\{|\operatorname{Cor}_{\text{Pearson}}(X^*,Y^*)|,\; | \operatorname{Cor}_{\text{Pearson}}((X^*)^2,(Y^*)^2)|\}} 
#' and 
#' \deqn{\rho_{L,1}(X,Y) = \max\{|\operatorname{Cor}_{\text{Pearson}}(X,Y)|,\; | \operatorname{Cor}_{\text{Pearson}}((X^*)^2,(Y^*)^2)|\}} 
#' the Lancaster correlation coefficient and the linear Lancaster correlation coefficient, respectively.
#' Two estimation methods are supported:
#'
#' \itemize{
#'   \item \strong{Linear estimator for \eqn{\bold{\rho_{L,1}}}} (\code{type = "linear"}): Consider \eqn{\rho_{L1} = \operatorname{Cor}_{\text{Pearson}}(X,Y)} and \eqn{\rho_{L2} = \operatorname{Cor}_{\text{Pearson}}((X^*)^2,(Y^*)^2)}.
#'   Let \eqn{\hat\rho_{L1}} be the sample Pearson correlation and \eqn{\hat\rho_{L2}} the empirical correlation of the squares of the empirically standardized observations, and set
#'   \eqn{\hat\rho_{L,1} = \max\{\,|\hat\rho_{L1}|,\;|\hat\rho_{L2}|\,\}}.
#'
#'   \item \strong{Rank-based estimator for \eqn{\bold{\rho_{L}}}} (\code{type = "rank"}): Consider \eqn{\rho_{R1} = \operatorname{Cor}_{\text{Pearson}}(X,Y)} and \eqn{\rho_{R2} = \operatorname{Cor}_{\text{Pearson}}((X^*)^2,(Y^*)^2)}.
#'   Let \eqn{Q_i} and \eqn{R_i} be the ranks of \eqn{X_i} and \eqn{Y_i}, within \eqn{X_1,...,X_n} or \eqn{Y_1,...,Y_n} respectively. Define
#'   \deqn{\hat\rho_{R1} = \frac{1}{n\,s_a^2}\sum_{i=1}^n a(Q_i)\,a(R_i),}
#'   \deqn{\hat\rho_{R2} = \frac{1}{n\,s_b^2}\sum_{i=1}^n \bigl(b(Q_i)-\bar b\bigr)\,\bigl(b(R_i)-\bar b\bigr),}
#'   where the scores are, for \eqn{j=1,...,n},
#'   \deqn{a(j) = \Phi^{-1}\!\Bigl(\frac{j}{n+1}\Bigr), \quad
#'   s_a^2 = \frac{1}{n}\sum_{j=1}^n\bigl(a(j)-\bar a\bigr)^2,} 
#'   \deqn{b(j)=a(j)^2,\quad
#'   \bar b=\frac{1}{n}\sum_{j=1}^n b(j), \quad
#'   s_b^2 = \frac{1}{n}\sum_{j=1}^n\bigl(b(j)-\bar b\bigr)^2.}
#'   Finally, the rank‐based Lancaster correlation is
#'   \deqn{\hat\rho_{L} = \max\bigl\{\,|\hat\rho_{R1}|, |\hat\rho_{R2}|\bigr\}.}
#' }
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
#' Holzmann, Klar (2024). "Lancester correlation - a new dependence measure linked to maximum correlation". \doi{https://doi.org/10.1111/sjos.12733}
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
