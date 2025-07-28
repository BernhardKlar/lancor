#' ACE permutation test of independence
#' 
#' @description 
#' Performs a permutation test of independence using ace in package acepack. Ace stands for alternating conditional expectations.
#'
#' @param x a numeric vector or a matrix-like object with two columns.
#' @param y NULL (default) or a vector with same length as x.
#' @param nperm number of permutations. 
#'  
#' @return A list containing the following components:
#' \item{ace}{the value of the test statistic}
#' \item{pval}{the p-value of the test}
#' 
#' @author Hajo Holzmann, Bernhard Klar
#' 
#' @references
#' Holzmann, Klar (2024). "Lancester correlation - a new dependence measure linked to maximum correlation". \url{https://arxiv.org/abs/2303.17872}
#' 
#' @seealso 
#' \code{\link{lcor.test}} for the Lancester correlation test.
#'
#'@examples 
#' n <- 200
#' x <- matrix(rnorm(n*2), n)
#' nu <- 2
#' y <- x / sqrt(rchisq(n, nu)/nu) #multivariate t
#' cor.test(y[,1], y[,2], method = "spearman")
#' ace.test(y)
#' 
#' @export
ace.test = function(x, y = NULL, nperm = 999) { #permutation test with ace
  if (is.data.frame(x)) x = as.matrix(x)
  if (!is.matrix(x) && is.null(y))
    stop("supply both 'x' and 'y' or a matrix-like 'x'")
  if (is.matrix(x)) {
    y = x[,2]
    x = x[,1]
  }
  a = acepack::ace(x, y)
  ace.cor =  as.vector( cor(a$tx, a$ty) )
  n = length(x)

  ts = ace.cor
  if (n <= 6) {
    nperm = arrangements::npermutations(n) #use all permutations
    perm = arrangements::permutations(x)
  } else {
    perm = arrangements::permutations(x, nsample = nperm)
  }
  tp = rep(0, nperm)
  for (i in 1:nperm) {
    a = acepack::ace(perm[i,], y)
    tp[i] = as.vector( cor(a$tx, a$ty) )
  }
  pval = (sum(tp > ts) + 1) / (nperm + 1)
  return(list(ace = ace.cor, pval = pval))
}
